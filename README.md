# Coordinates of Sights and Sounds

## Develop
```bash
npx parcel src/index.html
```

## Build
```bash
rm -r dist/
npx parcel build src/index.html
```

## Deploy
```bash
git add ./dist
git commit -m"vx.x.x"
git subtree push --prefix dist origin master
```
