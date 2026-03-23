#!/bin/bash

#echo -e '\n$ git pull origin'
#git pull origin main

echo -e '\n$ git add .'
git add .

#echo -e '\n$ git commit -m "Merge remote changes"'
#git commit -m "Merge remote changes"

echo -e '\n$ git commit -m "Update COMET"'
git commit -m "Update COMET"

echo -e '\n$ git push -u origin main'
git push -u origin main
