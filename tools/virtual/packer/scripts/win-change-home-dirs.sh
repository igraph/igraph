set -x

mkpasswd -l -p "$(cygpath $(cygpath -dH))" > /etc/passwd
