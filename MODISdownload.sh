#!/bin/bash

GREP_OPTIONS=''

cookiejar=$(mktemp cookies.XXXXXXXXXX)
netrc=$(mktemp netrc.XXXXXXXXXX)
chmod 0600 "$cookiejar" "$netrc"
function finish {
  rm -rf "$cookiejar" "$netrc"
}

trap finish EXIT
WGETRC="$wgetrc"

prompt_credentials() {
    echo "Enter your Earthdata Login or other provider supplied credentials"
    read -p "Username (luc.villandre): " username
    username=${username:-luc.villandre}
    read -s -p "Password: " password
    echo "machine urs.earthdata.nasa.gov login $username password $password" >> $netrc
    echo
}

exit_with_error() {
    echo
    echo "Unable to Retrieve Data"
    echo
    echo $1
    echo
    echo "https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.31/MOD11A1.A2012152.h25v06.006.2016112001610.hdf"
    echo
    exit 1
}

prompt_credentials
  detect_app_approval() {
    approved=`curl -s -b "$cookiejar" -c "$cookiejar" -L --max-redirs 2 --netrc-file "$netrc" https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.31/MOD11A1.A2012152.h25v06.006.2016112001610.hdf -w %{http_code} | tail  -1`
    if [ "$approved" -ne "302" ]; then
        # User didn't approve the app. Direct users to approve the app in URS
        exit_with_error "Please ensure that you have authorized the remote application by visiting the link below "
    fi
}

setup_auth_curl() {
    # Firstly, check if it require URS authentication
    status=$(curl -s -z "$(date)" -w %{http_code} https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.31/MOD11A1.A2012152.h25v06.006.2016112001610.hdf | tail -1)
    if [[ "$status" -ne "200" && "$status" -ne "304" ]]; then
        # URS authentication is required. Now further check if the application/remote service is approved.
        detect_app_approval
    fi
}

setup_auth_wget() {
    # The safest way to auth via curl is netrc. Note: there's no checking or feedback
    # if login is unsuccessful
    touch ~/.netrc
    chmod 0600 ~/.netrc
    credentials=$(grep 'machine urs.earthdata.nasa.gov' ~/.netrc)
    if [ -z "$credentials" ]; then
        cat "$netrc" >> ~/.netrc
    fi
}

fetch_urls() {
  if command -v curl >/dev/null 2>&1; then
      setup_auth_curl
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        curl -f -b "$cookiejar" -c "$cookiejar" -L --netrc-file "$netrc" -g -o $stripped_query_params -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
      done;
  elif command -v wget >/dev/null 2>&1; then
      # We can't use wget to poke provider server to get info whether or not URS was integrated without download at least one of the files.
      echo
      echo "WARNING: Can't find curl, use wget instead."
      echo "WARNING: Script may not correctly identify Earthdata Login integrations."
      echo
      setup_auth_wget
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        wget --load-cookies "$cookiejar" --save-cookies "$cookiejar" --output-document $stripped_query_params --keep-session-cookies -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
      done;
  else
      exit_with_error "Error: Could not find a command-line downloader.  Please install curl or wget"
  fi
}

fetch_urls <<'EDSCEOF'
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.31/MOD11A1.A2012152.h25v06.006.2016112001610.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.31/MOD11A1.A2012152.h25v07.006.2016112001603.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.31/MOD11A1.A2012152.h24v06.006.2016112001611.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.31/MOD11A1.A2012152.h24v07.006.2016112001605.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.30/MOD11A1.A2012151.h25v06.006.2016111234452.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.30/MOD11A1.A2012151.h24v06.006.2016111234426.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.30/MOD11A1.A2012151.h25v07.006.2016111234444.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.30/MOD11A1.A2012151.h24v07.006.2016111234426.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.29/MOD11A1.A2012150.h24v06.006.2016111221519.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.29/MOD11A1.A2012150.h25v07.006.2016111221529.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.29/MOD11A1.A2012150.h25v06.006.2016111221532.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.29/MOD11A1.A2012150.h24v07.006.2016111221515.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.28/MOD11A1.A2012149.h24v07.006.2016111214106.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.28/MOD11A1.A2012149.h25v07.006.2016111214101.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.28/MOD11A1.A2012149.h24v06.006.2016111214110.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.28/MOD11A1.A2012149.h25v06.006.2016111214111.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.27/MOD11A1.A2012148.h25v06.006.2016111200628.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.27/MOD11A1.A2012148.h25v07.006.2016111200614.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.27/MOD11A1.A2012148.h24v07.006.2016111200610.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.27/MOD11A1.A2012148.h24v06.006.2016111200619.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.26/MOD11A1.A2012147.h25v07.006.2016111193552.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.26/MOD11A1.A2012147.h25v06.006.2016111193546.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.26/MOD11A1.A2012147.h24v07.006.2016111193536.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.26/MOD11A1.A2012147.h24v06.006.2016111193550.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.25/MOD11A1.A2012146.h24v07.006.2016111172446.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.25/MOD11A1.A2012146.h24v06.006.2016111172459.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.25/MOD11A1.A2012146.h25v07.006.2016111172453.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.25/MOD11A1.A2012146.h25v06.006.2016111172507.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.24/MOD11A1.A2012145.h25v07.006.2016111124123.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.24/MOD11A1.A2012145.h24v07.006.2016111124111.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.24/MOD11A1.A2012145.h25v06.006.2016111124130.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.24/MOD11A1.A2012145.h24v06.006.2016111124119.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.23/MOD11A1.A2012144.h24v06.006.2016111111715.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.23/MOD11A1.A2012144.h25v06.006.2016111111709.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.23/MOD11A1.A2012144.h25v07.006.2016111111703.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.23/MOD11A1.A2012144.h24v07.006.2016111111701.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.22/MOD11A1.A2012143.h24v07.006.2016111110314.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.22/MOD11A1.A2012143.h25v07.006.2016111110320.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.22/MOD11A1.A2012143.h25v06.006.2016111110321.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.22/MOD11A1.A2012143.h24v06.006.2016111110326.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.21/MOD11A1.A2012142.h25v06.006.2016111095219.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.21/MOD11A1.A2012142.h24v06.006.2016111095215.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.21/MOD11A1.A2012142.h24v07.006.2016111095210.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.21/MOD11A1.A2012142.h25v07.006.2016111095218.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.20/MOD11A1.A2012141.h24v07.006.2016111093525.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.20/MOD11A1.A2012141.h24v06.006.2016111093532.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.20/MOD11A1.A2012141.h25v07.006.2016111093529.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.20/MOD11A1.A2012141.h25v06.006.2016111093531.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.19/MOD11A1.A2012140.h25v06.006.2016111082041.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.19/MOD11A1.A2012140.h24v06.006.2016111082026.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.19/MOD11A1.A2012140.h24v07.006.2016111082018.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.19/MOD11A1.A2012140.h25v07.006.2016111082034.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.18/MOD11A1.A2012139.h24v06.006.2016111080655.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.18/MOD11A1.A2012139.h25v06.006.2016111080713.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.18/MOD11A1.A2012139.h24v07.006.2016111080655.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.18/MOD11A1.A2012139.h25v07.006.2016111080701.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.17/MOD11A1.A2012138.h24v06.006.2016111064728.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.17/MOD11A1.A2012138.h25v06.006.2016111064723.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.17/MOD11A1.A2012138.h25v07.006.2016111064721.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.17/MOD11A1.A2012138.h24v07.006.2016111064719.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.16/MOD11A1.A2012137.h24v07.006.2016111062903.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.16/MOD11A1.A2012137.h24v06.006.2016111062835.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.16/MOD11A1.A2012137.h25v07.006.2016111062929.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.16/MOD11A1.A2012137.h25v06.006.2016111062937.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.15/MOD11A1.A2012136.h25v06.006.2016111050513.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.15/MOD11A1.A2012136.h24v07.006.2016111050457.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.15/MOD11A1.A2012136.h25v07.006.2016111050506.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.15/MOD11A1.A2012136.h24v06.006.2016111050501.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.14/MOD11A1.A2012135.h25v07.006.2016111045539.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.14/MOD11A1.A2012135.h24v07.006.2016111045538.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.14/MOD11A1.A2012135.h24v06.006.2016111045544.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.14/MOD11A1.A2012135.h25v06.006.2016111045541.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.13/MOD11A1.A2012134.h25v07.006.2016111031828.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.13/MOD11A1.A2012134.h25v06.006.2016111031835.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.13/MOD11A1.A2012134.h24v07.006.2016111031832.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.13/MOD11A1.A2012134.h24v06.006.2016111031823.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.12/MOD11A1.A2012133.h25v06.006.2016111030405.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.12/MOD11A1.A2012133.h24v07.006.2016111030350.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.12/MOD11A1.A2012133.h25v07.006.2016111030400.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.12/MOD11A1.A2012133.h24v06.006.2016111030351.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.11/MOD11A1.A2012132.h25v06.006.2016111012626.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.11/MOD11A1.A2012132.h24v07.006.2016111012611.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.11/MOD11A1.A2012132.h25v07.006.2016111012619.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.11/MOD11A1.A2012132.h24v06.006.2016111012616.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.10/MOD11A1.A2012131.h24v07.006.2016111011529.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.10/MOD11A1.A2012131.h24v06.006.2016111011537.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.10/MOD11A1.A2012131.h25v07.006.2016111011540.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.10/MOD11A1.A2012131.h25v06.006.2016111011546.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.09/MOD11A1.A2012130.h24v06.006.2016111000003.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.09/MOD11A1.A2012130.h24v07.006.2016110235936.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.09/MOD11A1.A2012130.h25v07.006.2016110235945.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.09/MOD11A1.A2012130.h25v06.006.2016110235943.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.08/MOD11A1.A2012129.h25v07.006.2016110232205.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.08/MOD11A1.A2012129.h24v07.006.2016110232206.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.08/MOD11A1.A2012129.h25v06.006.2016110232210.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.08/MOD11A1.A2012129.h24v06.006.2016110232214.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.07/MOD11A1.A2012128.h24v06.006.2016110221632.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.07/MOD11A1.A2012128.h25v07.006.2016110221656.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.07/MOD11A1.A2012128.h25v06.006.2016110221636.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.07/MOD11A1.A2012128.h24v07.006.2016110221621.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.06/MOD11A1.A2012127.h25v06.006.2016110214313.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.06/MOD11A1.A2012127.h25v07.006.2016110214323.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.06/MOD11A1.A2012127.h24v06.006.2016110214309.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.06/MOD11A1.A2012127.h24v07.006.2016110214310.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.05/MOD11A1.A2012126.h24v07.006.2016110203235.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.05/MOD11A1.A2012126.h25v07.006.2016110203236.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.05/MOD11A1.A2012126.h25v06.006.2016110203242.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.05/MOD11A1.A2012126.h24v06.006.2016110203237.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.04/MOD11A1.A2012125.h24v06.006.2016110200603.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.04/MOD11A1.A2012125.h25v06.006.2016110200655.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.04/MOD11A1.A2012125.h25v07.006.2016110200651.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.04/MOD11A1.A2012125.h24v07.006.2016110200559.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.03/MOD11A1.A2012124.h24v06.006.2016110185435.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.03/MOD11A1.A2012124.h25v06.006.2016110185714.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.03/MOD11A1.A2012124.h25v07.006.2016110185703.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.03/MOD11A1.A2012124.h24v07.006.2016110185430.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.02/MOD11A1.A2012123.h24v06.006.2016110182805.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.02/MOD11A1.A2012123.h25v07.006.2016110182805.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.02/MOD11A1.A2012123.h24v07.006.2016110182756.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.02/MOD11A1.A2012123.h25v06.006.2016110182817.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.01/MOD11A1.A2012122.h25v07.006.2016110165126.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.01/MOD11A1.A2012122.h24v06.006.2016110165112.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.01/MOD11A1.A2012122.h24v07.006.2016110165111.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLT/MOD11A1.006/2012.05.01/MOD11A1.A2012122.h25v06.006.2016110165134.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.31/MYD11A1.A2012152.h25v07.006.2016112043920.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.31/MYD11A1.A2012152.h25v06.006.2016112043920.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.31/MYD11A1.A2012152.h24v07.006.2016112043914.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.31/MYD11A1.A2012152.h24v06.006.2016112043917.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.30/MYD11A1.A2012151.h24v06.006.2016112042457.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.30/MYD11A1.A2012151.h25v06.006.2016112042514.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.30/MYD11A1.A2012151.h25v07.006.2016112042509.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.30/MYD11A1.A2012151.h24v07.006.2016112042447.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.29/MYD11A1.A2012150.h25v07.006.2016112025258.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.29/MYD11A1.A2012150.h24v06.006.2016112025308.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.29/MYD11A1.A2012150.h24v07.006.2016112025322.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.29/MYD11A1.A2012150.h25v06.006.2016112025309.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.28/MYD11A1.A2012149.h25v06.006.2016112023543.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.28/MYD11A1.A2012149.h24v06.006.2016112023533.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.28/MYD11A1.A2012149.h25v07.006.2016112023537.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.28/MYD11A1.A2012149.h24v07.006.2016112023529.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.27/MYD11A1.A2012148.h24v06.006.2016112011208.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.27/MYD11A1.A2012148.h25v06.006.2016112011219.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.27/MYD11A1.A2012148.h24v07.006.2016112011158.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.27/MYD11A1.A2012148.h25v07.006.2016112011216.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.26/MYD11A1.A2012147.h24v07.006.2016112004954.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.26/MYD11A1.A2012147.h24v06.006.2016112005002.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.26/MYD11A1.A2012147.h25v07.006.2016112004957.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.26/MYD11A1.A2012147.h25v06.006.2016112005007.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.25/MYD11A1.A2012146.h24v06.006.2016111232746.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.25/MYD11A1.A2012146.h25v07.006.2016111232735.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.25/MYD11A1.A2012146.h25v06.006.2016111232742.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.25/MYD11A1.A2012146.h24v07.006.2016111232733.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.24/MYD11A1.A2012145.h24v06.006.2016111230913.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.24/MYD11A1.A2012145.h25v06.006.2016111230832.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.24/MYD11A1.A2012145.h24v07.006.2016111230900.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.24/MYD11A1.A2012145.h25v07.006.2016111230822.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.23/MYD11A1.A2012144.h25v07.006.2016111213428.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.23/MYD11A1.A2012144.h24v07.006.2016111213425.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.23/MYD11A1.A2012144.h24v06.006.2016111213434.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.23/MYD11A1.A2012144.h25v06.006.2016111213439.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.22/MYD11A1.A2012143.h24v06.006.2016111212056.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.22/MYD11A1.A2012143.h24v07.006.2016111212045.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.22/MYD11A1.A2012143.h25v07.006.2016111212041.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.22/MYD11A1.A2012143.h25v06.006.2016111212055.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.21/MYD11A1.A2012142.h25v07.006.2016111193811.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.21/MYD11A1.A2012142.h25v06.006.2016111193812.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.21/MYD11A1.A2012142.h24v06.006.2016111193810.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.21/MYD11A1.A2012142.h24v07.006.2016111193809.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.20/MYD11A1.A2012141.h24v06.006.2016111193155.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.20/MYD11A1.A2012141.h25v06.006.2016111193158.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.20/MYD11A1.A2012141.h25v07.006.2016111193151.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.20/MYD11A1.A2012141.h24v07.006.2016111193151.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.19/MYD11A1.A2012140.h24v06.006.2016111122839.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.19/MYD11A1.A2012140.h25v06.006.2016111122842.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.19/MYD11A1.A2012140.h25v07.006.2016111122848.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.19/MYD11A1.A2012140.h24v07.006.2016111122826.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.18/MYD11A1.A2012139.h24v07.006.2016111122059.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.18/MYD11A1.A2012139.h24v06.006.2016111122104.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.18/MYD11A1.A2012139.h25v07.006.2016111122058.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.18/MYD11A1.A2012139.h25v06.006.2016111122103.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.17/MYD11A1.A2012138.h25v07.006.2016111072359.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.17/MYD11A1.A2012138.h24v07.006.2016111072340.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.17/MYD11A1.A2012138.h25v06.006.2016111072358.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.17/MYD11A1.A2012138.h24v06.006.2016111072349.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.16/MYD11A1.A2012137.h24v07.006.2016111061305.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.16/MYD11A1.A2012137.h25v07.006.2016111061322.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.16/MYD11A1.A2012137.h24v06.006.2016111061317.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.16/MYD11A1.A2012137.h25v06.006.2016111061326.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.15/MYD11A1.A2012136.h25v06.006.2016111053002.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.15/MYD11A1.A2012136.h24v06.006.2016111052957.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.15/MYD11A1.A2012136.h25v07.006.2016111052956.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.15/MYD11A1.A2012136.h24v07.006.2016111052945.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.14/MYD11A1.A2012135.h25v07.006.2016111042941.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.14/MYD11A1.A2012135.h25v06.006.2016111042942.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.14/MYD11A1.A2012135.h24v06.006.2016111042938.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.14/MYD11A1.A2012135.h24v07.006.2016111042925.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.13/MYD11A1.A2012134.h25v06.006.2016111035435.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.13/MYD11A1.A2012134.h25v07.006.2016111035425.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.13/MYD11A1.A2012134.h24v06.006.2016111035443.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.13/MYD11A1.A2012134.h24v07.006.2016111035425.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.12/MYD11A1.A2012133.h24v07.006.2016111030958.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.12/MYD11A1.A2012133.h25v06.006.2016111031003.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.12/MYD11A1.A2012133.h24v06.006.2016111030955.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.12/MYD11A1.A2012133.h25v07.006.2016111031000.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.11/MYD11A1.A2012132.h25v07.006.2016111021858.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.11/MYD11A1.A2012132.h24v06.006.2016111021836.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.11/MYD11A1.A2012132.h25v06.006.2016111021857.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.11/MYD11A1.A2012132.h24v07.006.2016111021848.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.10/MYD11A1.A2012131.h25v06.006.2016111010715.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.10/MYD11A1.A2012131.h25v07.006.2016111010710.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.10/MYD11A1.A2012131.h24v07.006.2016111010700.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.10/MYD11A1.A2012131.h24v06.006.2016111010709.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.09/MYD11A1.A2012130.h25v07.006.2016111001830.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.09/MYD11A1.A2012130.h25v06.006.2016111001839.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.09/MYD11A1.A2012130.h24v07.006.2016111001817.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.09/MYD11A1.A2012130.h24v06.006.2016111001826.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.08/MYD11A1.A2012129.h25v07.006.2016110231739.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.08/MYD11A1.A2012129.h24v06.006.2016110231737.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.08/MYD11A1.A2012129.h24v07.006.2016110231730.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.08/MYD11A1.A2012129.h25v06.006.2016110231746.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.07/MYD11A1.A2012128.h24v07.006.2016110224608.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.07/MYD11A1.A2012128.h25v06.006.2016110224624.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.07/MYD11A1.A2012128.h25v07.006.2016110224621.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.07/MYD11A1.A2012128.h24v06.006.2016110224619.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.06/MYD11A1.A2012127.h25v06.006.2016110213716.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.06/MYD11A1.A2012127.h24v06.006.2016110213709.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.06/MYD11A1.A2012127.h24v07.006.2016110213700.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.06/MYD11A1.A2012127.h25v07.006.2016110213712.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.05/MYD11A1.A2012126.h25v06.006.2016110210141.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.05/MYD11A1.A2012126.h25v07.006.2016110210135.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.05/MYD11A1.A2012126.h24v06.006.2016110210138.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.05/MYD11A1.A2012126.h24v07.006.2016110210122.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.04/MYD11A1.A2012125.h24v07.006.2016110194436.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.04/MYD11A1.A2012125.h25v06.006.2016110194454.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.04/MYD11A1.A2012125.h24v06.006.2016110194445.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.04/MYD11A1.A2012125.h25v07.006.2016110194456.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.03/MYD11A1.A2012124.h25v06.006.2016110193020.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.03/MYD11A1.A2012124.h25v07.006.2016110193028.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.03/MYD11A1.A2012124.h24v06.006.2016110193022.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.03/MYD11A1.A2012124.h24v07.006.2016110193013.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.02/MYD11A1.A2012123.h25v07.006.2016110175005.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.02/MYD11A1.A2012123.h24v06.006.2016110175015.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.02/MYD11A1.A2012123.h25v06.006.2016110175005.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.02/MYD11A1.A2012123.h24v07.006.2016110175007.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.01/MYD11A1.A2012122.h25v06.006.2016110174214.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.01/MYD11A1.A2012122.h24v06.006.2016110174114.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.01/MYD11A1.A2012122.h24v07.006.2016110174100.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Dal_E/MOLA/MYD11A1.006/2012.05.01/MYD11A1.A2012122.h25v07.006.2016110174210.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_C/MOTA/MCD12Q1.006/2012.01.01/MCD12Q1.A2012001.h24v07.006.2018146011412.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_C/MOTA/MCD12Q1.006/2012.01.01/MCD12Q1.A2012001.h25v06.006.2018146011453.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_C/MOTA/MCD12Q1.006/2012.01.01/MCD12Q1.A2012001.h25v07.006.2018146011505.hdf
https://e4ftl01.cr.usgs.gov//MODV6_Cmp_C/MOTA/MCD12Q1.006/2012.01.01/MCD12Q1.A2012001.h24v06.006.2018146011403.hdf
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N19E074.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N20E073.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N18E073.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N18E072.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N18E075.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N19E073.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N19E075.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N18E074.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N17E074.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N20E074.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N20E075.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N20E072.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N19E072.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N17E075.zip
https://e4ftl01.cr.usgs.gov//ASTER_B/ASTT/ASTGTM.003/2000.03.01/ASTGTMV003_N17E073.zip
EDSCEOF