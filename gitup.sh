#!/bin/bash

echo -e '\n$ git pull origin'
git pull origin main

echo -e '\n$ git add .'
git add .

echo -e '\n$ git commit -m "Merge remote changes"'
git commit -m "Merge remote changes"

echo -e '\n$ git push origin main'
git push origin main
