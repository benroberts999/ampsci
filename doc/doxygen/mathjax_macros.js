MathJax.Hub.Config({
  TeX: {
    Macros: {
      vb:       ["{\\boldsymbol{#1}}", 1],
      bra:      ["{\\langle #1 |}", 1],
      ket:      ["{| #1 \\rangle}", 1],
      matel:    ["{\\langle #1 | #2 | #3 \\rangle}", 3],
      redmatel: ["{\\langle #1 \\| #2 \\| #3 \\rangle}", 3],
      threej:   ["{\\begin{pmatrix}#1&#2&#3\\\\#4&#5&#6\\end{pmatrix}}", 6],
      sixj:     ["{\\begin{Bmatrix}#1&#2&#3\\\\#4&#5&#6\\end{Bmatrix}}", 6]
    }
  }
});
