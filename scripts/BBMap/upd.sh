#!/bin/sh

mkdir t tt 2>/dev/null

# extract BBMap tgz to ./t/ first.
tar -xzv -C ./t --strip-components=1 -f $1

rm -fr current && mv t/current .
rm -fr docs && mv t/docs .
rm -fr config && mv t/config .
rm -fr jni && mv t/jni .
rm -fr resources && mv t/resources .
rm -fr sh && mkdir sh && mv t/*.sh sh/
mv t/build.xml . && mv t/license.txt .
mv t/README.md tt/

sed -i.bak 's/^# *//g' tt/README.md
cat tt/README.md
ln -s ../current sh/
ls -la t tt sh/current

git add sh docs resources current config jni

echo git commit . -m "'Extract Version 3?.?? from `basename $1`'"
echo git tag -a v35.85 -F tt/README.md
# git push && git push --tags

