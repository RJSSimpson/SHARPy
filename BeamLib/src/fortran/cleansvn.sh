#!/bin/bash
for d in */*/*/*/ ; do 
  pushd $d
  sudo rm -r .svn
  popd
done
