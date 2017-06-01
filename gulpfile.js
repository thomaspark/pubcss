var gulp = require('gulp');
var browserSync = require('browser-sync').create();
var sass = require('gulp-sass');

// Static Server + watching scss/html files
gulp.task('serve', ['sass'], function() {

    browserSync.init({
        server: "./",
        directory: true
    });

    gulp.watch("./src/**/*.scss", ['sass']);
    // gulp.watch("app/*.html").on('change', browserSync.reload);
});

// Compile sass into CSS & auto-inject into browsers
gulp.task('sass', function() {
    return gulp.src("./src/*.scss")
        .pipe(sass())
        .pipe(gulp.dest("./dist/css/"))
        .pipe(browserSync.stream());
});

gulp.task('default', ['serve']);
