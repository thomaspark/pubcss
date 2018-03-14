const gulp = require('gulp');
const gutil = require('gulp-util');
// Developer
const browserSync = require('browser-sync').create();
// CSS Imports
const sass = require('gulp-sass');
const autoprefixer = require('gulp-autoprefixer');
const cleanCss = require('gulp-clean-css');
// PDFs
const puppeteerPdf = require('./gulp-puppeteer-pdf.js');

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

  //PDF
  pdf: {
    format: 'Letter',
    margin: {
      top: '1cm',
      right: '1cm',
      bottom: '1cm',
      left: '1cm'
    }
  },
  // BrowserSync
  browserSync: {
    server: {
        baseDir: "./dist/",
        directory: true
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
gulp.task('images', () => {
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
    .pipe(gulp.dest(`${setting.folder.dist}`));
});

// PDF
gulp.task('pdf', () => {
  gulp.src(`${setting.folder.dist}**/*.html`)
    .pipe(puppeteerPdf(setting.pdf));
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
  // Watch for updated files in the dist folder and trigger browserSync to
  // refresh.
  gulp.watch(`${setting.folder.dist}**/*.{css,html,png,gif,jpg,jpeg,svg,bmp}`, ['reload']);
  // Watch for changes in source folder and build file
  gulp.watch(`${setting.folder.src}**/*.scss`, ['css', 'pdf']);
  gulp.watch(`${setting.folder.src}**/*.html`, ['html', 'pdf']);
  gulp.watch(`${setting.folder.src}**/*.{png,gif,jpg,jpeg,svg,bmp}`, ['images', 'pdf']);
});

////////////////////////////////////////////////////////////////////////////////
// Meta Tasks
////////////////////////////////////////////////////////////////////////////////

// Build
gulp.task('build', ['html', 'images', 'css', 'pdf']);
// PDF must come last, becasue it's dependent on Final HTML/CSS actually
// exsisting.

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
    Builds documents, and generates PDF.

'npm run develop'
    Watches and builds

`);
});
