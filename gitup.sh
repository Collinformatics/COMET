#!/bin/bash

echo -e '\ngit pull origin'
git pull origin main

echo -e '\ngit add .'
git add .

echo -e '\ngit commit -m "Merge remote changes"'
git commit -m "Merge remote changes"

echo -e '\ngit push origin main'
git push origin main
