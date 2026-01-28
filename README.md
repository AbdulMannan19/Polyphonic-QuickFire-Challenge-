conda env create --prefix ./env -f environment.yml
conda activate ./env
conda env update --prefix ./env -f environment.yml --prune


git status
git add .
git commit -m "quick commit"
git push
