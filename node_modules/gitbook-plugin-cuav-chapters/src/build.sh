# Cleanup folder
rm -rf book

# Recreate folder
mkdir book

# Compile JS
uglifyjs src/js/cuav-chapters.js -o book/cuav-chapters.js -mc

# Compile Website CSS
node-sass --output-style compressed src/scss/cuav-chapters.scss > book/cuav-chapters.css
