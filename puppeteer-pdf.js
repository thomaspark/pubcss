'use strict';
const gutil = require('gulp-util');
const through = require('through2');
const puppeteer = require('puppeteer');

module.exports = (options) => through.obj((file, enc, cb) => {
  console.log(options);
  const makePdf = async (file) => {
    const pdfPath = gutil.replaceExtension(file.path, '.pdf');
    const browser = await puppeteer.launch();
    const page = await browser.newPage();
    await page.goto(`file://${file.path}`, {waitUntil: 'networkidle2'});
    await page.pdf({path: pdfPath, format: options.format, margin: options.margin});
    await browser.close();
  }

  if (file.isNull()) {
    cb(null, file);
    return;
  }

  if (file.isStream()) {
    cb(new gutil.PluginError('PuppeteerPdf™', 'Streaming not supported'));
    return;
  }
  makePdf(file);

  cb(null, file);
});
