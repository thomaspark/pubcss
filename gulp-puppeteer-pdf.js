'use strict';
const gutil = require('gulp-util');
const through = require('through2');
const puppeteer = require('puppeteer');

module.exports = (options) => through.obj((file, enc, cb) => {
  //A general function to save the PDF.
  const makePdf = async (file) => {
    // Made PDF paths.
    const pdfPath = gutil.replaceExtension(file.path, '.pdf');
    // Turn on Chrome
    const browser = await puppeteer.launch();
    // Make a new page!
    const page = await browser.newPage();
    // In the new page, go to the file.
    await page.goto(`file://${file.path}`, {waitUntil: 'networkidle2'});
    //PDF the crap out of that page.
    await page.pdf({
      // Grab it from the File System
      path: pdfPath,
      // "Paper Size"
      format: options.format,
      // We're using margin here so we don't print to the to the edge.
      margin: options.margin
    });
    // Close the session.
    await browser.close();
}

  if (file.isNull()) {
    cb(null, file);
    return;
  }

  // No streams for you!
  if (file.isStream()) {
    cb(new gutil.PluginError('PuppeteerPdf™', 'Streaming not supported'));
    return;
  }
  makePdf(file);

  cb(null, file);
});
