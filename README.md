### An amazing project
# PubCSS


PubCSS is a library of CSS stylesheets and HTML templates for formatting academic publications for print and the web.

Currently, the following formats are supported:

* [ACM SIG Proceedings](http://www.acm.org/sigs/publications/proceedings-templates)
* [ACM SIGCHI Conference Papers](http://www.sigchi.org/publications/chipubform)
* [ACM SIGCHI Extended Abstracts](http://www.sigchi.org/publications/chipubform)
* [IEEE Conference Proceedings](http://www.ieee.org/conferences_events/conferences/publishing/templates.html)

Check out sample output [here](http://thomaspark.me/2015/01/pubcss-formatting-academic-publications-in-html-css/).


## Quick Start

1. Create an HTML file, or modify the relevant template, to add your content
2. Link your HTML to the relevant `pub.css` stylesheet
3. Build to PDF using Prince with the command `prince input.html output.pdf`

The only dependency is [Prince](http://www.princexml.com/), which is free for non-commercial use.

## Reference

### Sections

Sections and subsections are automatically numbered by following this pattern.


```html
<h1>Section Header</h1>
<section>
  <h2>Subsection Header</h2>
  <p>Lorem ipsum</p>
</section>
```

### Tables and Figures

Figures and tables are also numbered if you include a caption.

```html
<table>
  <tr><td>1</td><td>2</td></tr>
  <caption>Example Table</caption>
<table>

<figure>
  <img src="graph.png">
  <figcaption>Example Figure</figcaption>
</figure>
```

### References and Citations

References are also numbered. Be sure to assign them unique IDs.

```html
<cite id="nicole">Nicole. 2015. Title of paper. <em>Journal</em>, 4(3), 1-10.</cite>
```

Citations to the references make use of these IDs.

```html
<a href="#nicole"></a>
```

Multiple citations can be made in one set of brackets.

```html
<span class="cites"><a href="#beeker"></a><a href="#jackie"></a><a href="#kiwi"></a></span>
```

Sections, tables, and figures can also be referenced by adding a class.

```html
<a href="#intro" class="section"></a>
<a href="#example-table" class="table"></a>
<a href="#example-figure" class="figure"></a>
```

### Equations

Equations are also numbered. MathML is well-supported by Prince. For the web, you’ll need [MathJax](http://www.mathjax.org/) to render MathML properly in Chrome and Internet Explorer.

```html
<div class="equation">
  <math xmlns="http://www.w3.org/1998/Math/MathML">
    <mi>E</mi>
    <mo>=</mo>
    <mi>m</mi>
    <msup>
      <mi>c</mi>
      <mn>2</mn>
    </msup>
  </math>
</div>
```

### Footnotes

Footnotes are made within the body text, and are automatically moved to the bottom of the current page.

```html
<p>This is text.<span class="footnote">And this is a footnote.</span></p>
```

### Smart Quotes

Smart quotes can be used in lieu of straight quotes by enclosing the text like so. You can nest quotes within quotes.

```html
<q>To be or not be.</q>
```

### Utility Classes

Utility classes are also available to modify layout and counter behavior.

* `col-2`: divide the element into 2 columns
* `col-3`: divide the element into 3 columns
* `col-4`: divide the element into 4 columns
* `col-span`: span all of the columns of parent
* `col-break-before`: force column break before element
* `col-break-after`: force column break after element
* `page-break-before`: force page break before element
* `page-break-after`: force page break after element
* `counter-skip`: do not count this header
* `counter-reset`: reset counter to 0

## Customization

One of the major advantages of PubCSS is that you can use CSS to customize the style. All of the usual rules apply.

To create a new theme, you’ll want to dig into the LESS source. The most common changes can be made through `variables.less`, such as toggling page numbers and setting counter styles. The rest can be included in `custom.less`.
