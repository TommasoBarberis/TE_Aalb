#!/bin/bash

DATADIR=$1
SOCKETDIR=$2

: "${DATADIR:?Must supply the path to the mysql data directory as the first argument to this script}"
: "${SOCKETDIR:?Must supply the path to the mysql socket directory as the second argument to this script}"

echo "Using \"$DATADIR\" as writable mysql database directory"
echo "Will create mysql communication socket in \"$SOCKETDIR\""
echo
echo "If these directories do not exist they will be created. Otherwise the contents of existing directories will be used."

case $DATADIR in
    /*) DATADIR_PATH="$DATADIR" ;;
    *) DATADIR_PATH="$PWD/$DATADIR" ;; 
esac

case $SOCKETDIR in
    /*) SOCKETDIR_PATH="$SOCKETDIR" ;;
    *) SOCKETDIR_PATH="$PWD/$SOCKETDIR" ;; 
esac

mkdir -p $DATADIR_PATH
mkdir -p $SOCKETDIR_PATH


singularity instance start --bind ${DATADIR_PATH}:/var/lib/mysql --bind ${SOCKETDIR_PATH}:/run/mysqld PASTEC.sif mysql 
singularity run instance://mysql 