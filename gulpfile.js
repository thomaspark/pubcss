const gulp = require('gulp');
const gutil = require('gulp-util');
// Developer
const browserSync = require('browser-sync').create();
// CSS Imports
const sass = require('gulp-sass');
const autoprefixer = require('gulp-autoprefixer');
const cleanCss = require('gulp-clean-css');

////////////////////////////////////////////////////////////////////////////////
// Settings
////////////////////////////////////////////////////////////////////////////////

// Browser Versions to Transpile/Prefix to.
const browserVersions = ['ie >= 11', 'last 2 versions'];

// Build Settings
const setting = {

  // Default Directories
  folder: {
    src: 'src/',
    dist: 'dist/'
  },

  // Autoprefixer
  autoprefixer: {
    browsers: browserVersions,
    cascade: true
  },

  // BrowserSync
  browserSync: {
    server: {
        baseDir: "./dist/",
        directory: false
    },
    ghostMode: {
      clicks: true,
      forms: true,
      scrolls: true
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Builds
////////////////////////////////////////////////////////////////////////////////

// Html
gulp.task('html', () => {
  gulp.src(`${setting.folder.src}**/*.html`)
    // Write to Dist Folder
    .pipe(gulp.dest(setting.folder.dist));
});

// Images
gulp.task('html', () => {
  gulp.src(`${setting.folder.src}**/*.{png,gif,jpg,jpeg,svg,bmp}`)
    // Write to Dist Folder
    .pipe(gulp.dest(setting.folder.dist));
});

// CSS
gulp.task('css', () => {
  gulp.src(`${setting.folder.src}**/*.scss`)
    // Transpile Sass to CSS
    .pipe(sass().on('error', sass.logError))
    // Autoprefix
    .pipe(autoprefixer(setting.autoprefixer))
    // Minify CSS
    .pipe(cleanCss())
    // Writing to Dist folder
    .pipe(gulp.dest(`${setting.folder.dist}css/`));
});

////////////////////////////////////////////////////////////////////////////////
// Development
////////////////////////////////////////////////////////////////////////////////

// Serve
gulp.task('serve', () => {
    browserSync.init(setting.browserSync);
});

// Reload
gulp.task('reload', () => {
  browserSync.reload();
});

// Watch
gulp.task('watch', () => {
  // Watch for updated files in the dist folder and reload
  gulp.watch(`${setting.folder.dist}**/*.*`, ['reload']);
  // Watch for changes in source folder and build file
  gulp.watch(`${setting.folder.src}**/*.scss`, ['css']);
  gulp.watch(`${setting.folder.src}**/*.html`, ['html']);
});

////////////////////////////////////////////////////////////////////////////////
// Meta Tasks
////////////////////////////////////////////////////////////////////////////////

// Build
gulp.task('build', ['html', 'css']);

// Develop
gulp.task('develop', ['build', 'serve', 'watch']);

////////////////////////////////////////////////////////////////////////////////
// Default
////////////////////////////////////////////////////////////////////////////////

gulp.task('default', () => {
  console.log(
`

Build should be run with 'npm run [command]'
--------------------------------------------
Options:

'npm run build'
    Builds project with no developer features

'npm run develop'
    Watches and builds

`);
});
