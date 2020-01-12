#!/bin/bash -x

# For python 3.6 on ubuntu 16.04
# https://askubuntu.com/questions/865554/how-do-i-install-python-3-6-using-apt-get
if ! which python3.6
then
  sudo add-apt-repository --yes ppa:jonathonf/python-3.6
  sudo apt-get update
fi
sudo apt-get install --yes python3.6 python3.6-dev

# pip3.6 for installing packages specifically for python 3.6
# https://askubuntu.com/questions/889535/how-to-install-pip-for-python-3-6-on-ubuntu-16-10
if ! which pip3.6
then
  curl https://bootstrap.pypa.io/get-pip.py | python3.6 - --user
fi

# Install tre regexp-up-to-edit-distance library
sudo apt install agrep libtre-dev
git clone https://github.com/ahomansikka/tre
(
  cd tre/python3
  python3.6 setup.py install --user
)

pip3.6 install --user biopython
pip3.6 install --user tqdm
pip3.6 install --user zstandard
