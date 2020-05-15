# Install script

## Installing on you Ubuntu
In your ubuntu, from this directory, run:

```
./ALL.sh
```

If that fails, please open an issue (with full output) at https://github.com/anpc/sines-in-aging/issues

## Testing the install script, without touching your machine

As we use different versions of ubuntu (some via Windows Subsystem for Linux), it became hard to write an install script that will work for everyone.

We can use "containers" to test it reproducibly, over and over again, on various ubuntu versions.
(We don't currently recommend doing actual research in containers, only _testing_ the install script. The final goal is to install on your machine.)

### Installing docker (popular tool for working with containers)

```
docker run hello-world
```
If the above command worked, you already have a working Docker.  If not:
 
```
sudo apt install docker.io
sudo adduser $USER docker
```
Now logout and login again.
Run `groups` and make sure it shows `docker` as one of your groups.
Retry the above "hello-world" command.

### Running the tests
From this directory, run:
```
docker build --file=Dockerfile.on-ubuntu-16.04 .
docker build --file=Dockerfile.on-ubuntu-18.04 .
```

### Cleaning disk space

Docker saves filesystem snapshots from _all runs_ and will happily eat all your disk!  To clean up:
```
docker system prune
```