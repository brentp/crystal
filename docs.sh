(git checkout gh-pages && git pull origin master) || exit
git checkout master

set -e
cd doc && make html
cd -
ghp-import -np doc/build/html/
