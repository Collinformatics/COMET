#!/bin/bash

# Update a github with the terminal

echo -e '\n$ git pull origin'
git pull origin main

echo -e '\n$ git add .'
git add .

echo -e '\n$ git commit -m "Merge remote changes"'
git commit -m "Merge remote changes"

echo -e '\n$ git push origin main'
git push origin main

echo -e '\n$ git checkout update'
git checkout update

echo -e '\n$ git pull origin main --rebase'
git pull origin main --rebase

echo -e '\n$ git push origin update'
git push origin update
