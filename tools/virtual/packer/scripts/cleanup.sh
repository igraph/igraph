#!/bin/bash -x
echo "Cleaning up dhcp leases..."
rm /var/lib/dhcp/*

echo "Cleaning up udev rules..."
rm /etc/udev/rules.d/70-persistent-net.rules
mkdir /etc/udev/rules.d/70-persistent-net.rules
rm -rf /dev/.udev/
rm /lib/udev/rules.d/75-persistent-net-generator.rules

#apt cleanup
echo "Running apt-get remove kernel headers..."
apt-get -y remove linux-headers-$(uname -r)
echo "Running remove older kernel headers"
dpkg -l 'linux-*' | sed '/^ii/!d;/'"$(uname -r | sed "s/\(.*\)-\([^0-9]\+\)/\1/")"'/d;s/^[^ ]* [^ ]* \([^ ]*\).*/\1/;/[0-9]/!d' | xargs sudo apt-get -y purge
echo "Running apt-get clean..."
apt-get -y clean
echo "Running apt-get autoclean..."
apt-get -y autoclean
echo "Running apt-get remove..."
apt-get -y remove
echo "Running apt-get auto-remove..."
apt-get -y autoremove

echo "pre-up sleep 2" >> /etc/network/interfaces

#zero out disk space. Replacing free space with 0s makes the drive more easily compressed
echo "Zeroing out disk..."
dd if=/dev/zero of=/EMPTY bs=1M || true
rm -f /EMPTY
exit

