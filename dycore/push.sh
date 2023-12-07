#!bin/sh

echo "# Dycore with passive tracer and latent heat" > README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote remove origin
git remote add origin https://github.com/kuiper2000/Dycore.git
git push -u origin main


git init
git fetch origin master
git merge origin master
git remote add origin https://github.com/kuiper2000/website-hugo.git
git add .
git commit -m "Initial commit"
git push -u origin master --force
