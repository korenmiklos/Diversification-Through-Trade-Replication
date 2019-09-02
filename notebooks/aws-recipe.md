1. Go to [AWS console](https://eu-central-1.console.aws.amazon.com/ec2/v2/home?region=eu-central-1#LaunchTemplates:sort=launchTemplateId|)
2. Launch instance from template, Version 7.
3. Once [instance is running](https://eu-central-1.console.aws.amazon.com/ec2/v2/home?region=eu-central-1#Instances:instanceId=i-0892e3c38c5d6f7d7;sort=instanceId), copy its domain name to the clipboard. 
4. From the shell within GD folder, `ssh -i awsjulia.pem ubuntu@ec2-3-122-225-122.eu-central-1.compute.amazonaws.com`
5. Mount the SBS block:
```
lsblk
sudo mkdir /data
sudo mount /dev/nvme1n1 /data
```
6. Create a symlink to Julia and to private keys:
```
sudo ln -s /home/ubuntu/JuliaPro-0.6.2.1/julia /usr/local/bin/julia
ln -s /data/id_rsa /home/ubuntu/.ssh/id_rsa
```
7. Pull the most recent version of the repo:
```
cd /data/impvol-julia
git stash
git pull
```
8. Install missing components.
```
sudo apt install htop
make install
```
9. Run make: `make tables`

Or just run all these from a shell script

`bash <(curl -s https://raw.githubusercontent.com/korenmiklos/impvol-julia/master/notebooks/aws-julia.sh)`