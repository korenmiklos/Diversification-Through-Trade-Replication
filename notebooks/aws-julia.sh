screen
lsblk
# set up environment
sudo mkdir /data
sudo mount /dev/nvme1n1 /data
sudo ln -s /home/ubuntu/JuliaPro-0.6.2.1/julia /usr/local/bin/julia
ln -s /data/id_rsa /home/ubuntu/.ssh/id_rsa
cd /data/impvol-julia
# pull most recent version
git stash
yes | git pull
# install missing components
sudo apt install htop
make install
