#!/usr/bin/env node
/* Builds heatmap/heatmap-standalone.html — a single self-contained file with
   the stylesheet, data bundle, and app script all inlined. Zero dependencies:
   double-click it anywhere and it works. Useful for sharing the heat map.

   Run after tools/standardize.js (it consumes the generated data.bundle.js).
   Usage:  node tools/build-standalone.js [heatmap-dir] */

'use strict';
const fs = require('fs');
const path = require('path');

const HM = process.argv[2] || 'heatmap';
const read = f => fs.readFileSync(path.join(HM, f), 'utf8');

/* `</script` inside an inlined script would close the tag early — neutralize it.
   The backslash is inert inside JS strings (the only place it can occur here). */
const safe = js => js.replace(/<\/script/gi, '<\\/script');

/* Use function replacements: a string replacement would interpret `$&`, `$\``
   etc. inside the inlined code as substitution patterns and corrupt it. */
let html = read('index.html');
html = html.replace('<link rel="stylesheet" href="styles.css">',
  () => '<style>\n' + read('styles.css') + '\n</style>');
html = html.replace('<script src="data.bundle.js"></script>',
  () => '<script>\n' + safe(read('data.bundle.js')) + '\n</script>');
html = html.replace('<script src="app.js"></script>',
  () => '<script>\n' + safe(read('app.js')) + '\n</script>');

const OUT = path.join(HM, 'heatmap-standalone.html');
fs.writeFileSync(OUT, html);
console.log('Wrote ' + OUT + ' (' + Math.round(html.length / 1024) + ' KB, self-contained)');
