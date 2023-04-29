#!/bin/zsh
set -x
if [ "$1" != "up" ]&& [ "$1" != "down" ]
then
	echo "up or down"
	exit 1
fi

case "$2" in
	"ccom")
	server=zzirui@ccom-boom-login.ucsd.edu
	;;

	"tscc")
	server=zzirui@tscc-login.sdsc.edu
	;;

	"181")
	server=zzirui@math181.ucsd.edu
	remote=c25/
	;;

	"177")
	server=zzirui@math177.ucsd.edu
	remote=c25/
	;;

	"181admin")
	server=zziruiadmin@math181.ucsd.edu
	remote=c25/
	;;

	"177admin")
	server=zziruiadmin@math177.ucsd.edu
	remote=c25/
	;;

	*)
	echo "unknown server"
	exit 1
	;;
esac


# declare -A newmap
# newmap[/Users/zzirui/c25/]=/home/zzirui/c25/
cdir=$PWD/ #need trailing / for rsync

# for i in "${(k)newmap[@]}" #loop through map
# do
# 	if [[ "$cdir" =~ "$i" ]] #if path matches
# 	then
# 		remote=${cdir/$i/${newmap[$i]}} 
# 		successfind = true
# 	fi
# done

remote=$server:$remote #create remote dir

if [ "$1" = "up" ]
then
	dest=$remote
	source=$cdir
	filelist=(--exclude={'*.o','*.out','.git/*','results/*','meeting/*','doc/*','bvismtests/*'})
else
	source=$remote
	dest=$cdir
	filelist=(--include={'results/*','bvismtests/*'})
fi

# https://superuser.com/questions/360966/how-do-i-use-a-bash-variable-string-containing-quotes-in-a-command


extraopt=(${@:3})


echo "source = $source, dest = $dest"

# echo $filelist
rsync -rtuv ${source} ${dest} ${filelist[@]} ${extraopt[@]}