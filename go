set -x
./ADSB_Encoder.py  0xABCDEF 12.34 56.78 9999.0
dd if=Samples.iq8s of=Samples_256K.iq8s bs=4k seek=63
sudo hackrf_transfer -t Samples_256K.iq8s -f 1090000000 -s 2000000 -x 10
