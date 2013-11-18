#!/bin/bash
# Compiles and installs all the supported Python distributions using pyenv

VERSIONS="2.6.9 2.7.5 3.1.5 3.2.5 3.3.2"

cd ~
if [ ! -d python ]; then
    mkdir python
fi

if [ ! -d pyenv ]; then
    git clone git://github.com/yyuu/pyenv.git pyenv
else
    cd pyenv
    git pull
    cd ..
fi

echo 'export PYENV_ROOT="$HOME/pyenv"' > ~/pyenv/env.sh
echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/pyenv/env.sh
echo 'eval "$(pyenv init -)"' >> ~/pyenv/env.sh

if grep -q "pyenv/env.sh" ~/.bashrc; then
    # Nothing to do
    true
else
    echo "" >>~/.bashrc
    echo "# Initialize PyEnv" >>~/.bashrc
    echo "source pyenv/env.sh" >>~/.bashrc
fi

source pyenv/env.sh

for VER in ${VERSIONS}; do
    if pyenv versions | grep -F -q "${VER}" ; then
        # Version already installed
        true
    else
        pyenv install ${VER}
    fi
    BASEVER="`echo ${VER} | cut -d '.' -f 1-2`"
    rm -f "python/${BASEVER}"
    ln -s "../pyenv/versions/${VER}" "python/${BASEVER}"
done
pyenv rehash
