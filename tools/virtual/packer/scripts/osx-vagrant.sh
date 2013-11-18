#!/bin/sh
date > /etc/vagrant_box_build_time
OSX_VERS=$(sw_vers -productVersion | awk -F "." '{print $2}')

# Set computer/hostname
COMPNAME=vagrant-osx-10-${OSX_VERS}
scutil --set ComputerName ${COMPNAME}
scutil --set HostName ${COMPNAME}.vagrantup.com

# Installing vagrant keys
mkdir /Users/vagrant/.ssh
chmod 700 /Users/vagrant/.ssh
curl -k 'https://raw.github.com/mitchellh/vagrant/master/keys/vagrant.pub' > /Users/vagrant/.ssh/authorized_keys
chmod 600 /Users/vagrant/.ssh/authorized_keys
chown -R vagrant /Users/vagrant/.ssh
