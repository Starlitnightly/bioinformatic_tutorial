/*
	#
*/
var exec        = require('child_process').exec;
var gulp 		= require("gulp");
var concat      = require("gulp-concat");
var merge       = require("merge-stream");
var rename 		= require("gulp-rename");
var sass 		= require("gulp-sass");
var browserify  = require('gulp-browserify');
var bsync 		= require("browser-sync");
var path        = require('path');

/*
    ## Copy over files
*/

gulp.task( "copy", function(done) {

    // Copy over stuff

    var fonts = gulp.src( ["node_modules/font-awesome/fonts/**/"] )
    .pipe( gulp.dest("_assets/website/fonts/fontawesome/") )
    ;

    var img = gulp.dest("_assets/website/images/")
    ;

    var merged = merge( fonts, img );
    // merged.add(js);


    // return merged;
    done();

});


/*
    ## Build gitbook core JS
*/

gulp.task('core-js', function()
{
    return  gulp.src('src/js/core/index.js')
            .pipe(browserify({
                insertGlobals : true
            }))
            .pipe(rename('gitbook.js'))
            .pipe(gulp.dest('./_assets/website/'))
});

/*
    ## Build theme JS
*/

gulp.task('theme-js', function()
{
    return  gulp.src('src/js/theme/index.js')
            .pipe(browserify({
                insertGlobals : true
            }))
            .pipe(rename('theme.js'))
            .pipe(gulp.dest('./_assets/website/'))
});


/*
    ## Compile the SASS into a single CSS file
*/

gulp.task('sass', function() {

    console.log( "sassing" );
    // Compile
    return  gulp.src('./src/scss/website.scss')
            .pipe( sass().on( 'error', sass.logError ) )
            .pipe(rename('style.css'))
            .pipe( gulp.dest( './_assets/website/' ) )
    ;
});

gulp.task("reload", function(done)
{
    console.log("Reloading gitbook....")
    exec('cd ../starter-kit/documentation/ && gitbook build', function (err, stdout, stderr)
    {
        console.log(stdout);
        console.log(stderr);
        bsync.reload();
    });

    done();

});

/*
    ## Development web server and file watcher
*/

gulp.task("server", function(done) {

    bsync.init(
    {
        server: '../starter-kit/documentation/_book',
        open: false
    });

    // Watch files
    gulp.watch( 'src/js/core/**/*.js', gulp.series('core-js', 'reload') );
    gulp.watch( 'src/js/theme/**/*.js', gulp.series('theme-js', 'reload') );
    gulp.watch( 'src/**/*.scss', gulp.series('sass', 'reload') );
    gulp.watch( '_layouts/**/*.html', gulp.series('reload') );

    done();

});



/*
    # Define main gulp tasks
*/

gulp.task("build", gulp.series( "copy", "core-js", "theme-js", "sass" ) );
gulp.task("default", gulp.series( "build", "server", "reload" ) );
