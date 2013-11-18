#set -x

VAGRANT_HOME=/cygdrive/c/Users/vagrant
# Install ssh certificates
mkdir $VAGRANT_HOME/.ssh
chmod 700 $VAGRANT_HOME/.ssh
cd $VAGRANT_HOME/.ssh
wget --no-check-certificate 'https://raw.github.com/mitchellh/vagrant/master/keys/vagrant.pub' -O authorized_keys
chown -R vagrant $VAGRANT_HOME/.ssh
cd ..

cd $VAGRANT_HOME

ZIP_INSTALL=7z922-x64.msi
VBOX_INSTALL=VBoxWindowsAdditions-amd64.exe
VMWARE_INSTALL=setup64.exe

# 7zip will allow us to extract a file from an ISO
wget http://downloads.sourceforge.net/sevenzip/$ZIP_INSTALL
msiexec /qb /i $ZIP_INSTALL

if [ -f VBoxGuestAdditions.iso ]; then
    # Extract the installer from the ISO
    /cygdrive/c/Program\ Files/7-Zip/7z.exe x VBoxGuestAdditions.iso $VBOX_INSTALL

    # Mark Oracle as a trusted installer
    certutil -addstore -f "TrustedPublisher" $(cygpath -d /cygdrive/a/oracle-cert.cer)

    # Install the Virtualbox Additions
    ./$VBOX_INSTALL /S

    # Cleanup
    rm -f VBoxGuestAdditions.iso
    rm -f $VBOX_INSTALL.exe
elif [ -f windows.iso ]; then
    # Extract the installer from the ISO
    /cygdrive/c/Program\ Files/7-Zip/7z.exe x windows.iso $VMWARE_INSTALL

    # Install VMware tools
    ./$VMWARE_INSTALL /S /v "/qn REBOOT=R ADDLOCAL=ALL" || true

    # Cleanup
    rm -f windows.iso
    rm -f $VMWARE_INSTALL
fi

# Cleanup
rm -f $ZIP_INSTALL
