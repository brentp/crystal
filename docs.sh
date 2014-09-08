set -e
cd doc && make html
cd -
ghp-import -np doc/build/html/
