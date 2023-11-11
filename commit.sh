rm *.user

cd Main
rm knitro.log
rm imgui.ini
rm *.trajopt
rm *.gif
rm *.mp2
cd .. 

git add --all
git commit -m "$1"
git push -u origin $2
