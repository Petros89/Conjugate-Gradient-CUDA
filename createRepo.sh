#!/bin/sh

# This script initiates a new repository name on Github

repo_name=$1

test -z $repo_name && echo "Repo name required." 1>&2 && exit 1

curl -u "Petros89" https://api.github.com/user/repos -d "{\"name\":\"$repo_name\"}"


# This script explains the steps to upload the files to the new directory on Github

# 1) Go inside the directory that the files you want to upload reside (e.g cd ~/myfolder/)

# 2) git init

# 3) git add . && git commit -m "your commit message"

# 4) git remote add origin "copy the Link of the repository on Github"

# 5) git push -u origin master -f


###########
#DEBUGGING#
###########

# If you face problems do the followings for username and e-mail veification

# git config --global user.email "your github email"
# git config --global user.name "your github username"
