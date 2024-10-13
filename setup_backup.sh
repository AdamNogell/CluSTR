#!/bin/bash

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

brew install bwa
brew install samtools
brew install gcc

brew install pipx
pipx install cutadapt
source ~/.bashrc

git clone https://github.com/weizhongli/cdhit.git
cd cdhit || return
make CC=/opt/homebrew/Cellar/gcc/14.2.0/bin/g++-14
sudo mv cd-hit /usr/local/bin
sudo mv cd-hit-est /usr/local/bin
cd .. || return
rm -r cdhit

git clone https://github.com/Debian/fastx-toolkit.git

brew install python
echo 'alias python=python3' >> ~/.bash_profile
echo 'alias pip=pip3' >> ~/.bash_profile
source ~/.bash_profile

pip install biopython