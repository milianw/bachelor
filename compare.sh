make || exit
echo "mw:"
time ./hs_mw &>/dev/null
echo "hossam:"
time ./hs_hossam &>/dev/null
