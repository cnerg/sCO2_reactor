#!/bin/bash

filename=$1

echo $filename
# Setup
mkdir aaswenson

# copy input to tmp directory
cp $filename aaswenson

cd aaswenson

# copy mcnp/data here
cp /mnt/gluster/aaswenson/test.tar.gz .
tar -xf test.tar.gz
rm test.tar.gz


# run
./mcnp6 i= $filename o= ../$filename.o x= testxsdir

cd 
rm -rf aaswenson
