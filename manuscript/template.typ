#let project(title: "", abstract: [], authors: (), body) = {
  set document(author: authors.map(a => a.name), title: title)
  set page(numbering: "1", number-align: start)
  set text(font: "Linux Libertine", lang: "en")

  // Title row.
  align(left)[
    #block(text(weight: 700, 2em, title)) \
  ]
  
  grid(
    ..authors.map(author => text()[#author.name,]),
  )

  // Abstract.
  pad(
    x: 1em,
    top: 1em,
    bottom: 1.1em,
    align(left)[
      #heading(
        outlined: false,
        numbering: none,
        text(0.8em, smallcaps[Abstract]),
      )
      #abstract
    ],
  )


  body
}

#let showtitle(title) = {    
    text(size: 26pt, fill: rgb("#222"), weight: 600)[#par(text(title), leading: 0.3em)]  
}

#let showabstract(abstract) = {
  heading(
        outlined: false,
        numbering: none,
        text(
          0.75em,     
          fill: rgb("#555"),
          [Abstract]
        ),
  )
  par(first-line-indent: 0em)[
    #text(
      size: 0.95em
    )[#abstract]
  ]
}


#let showauthors(authors, institutions, corresponding_email) = {

for person in authors.keys() [
  #text(size: 1em)[#person#for x in authors.at(person).affiliations {  
  let pos = institutions.keys().position(y => y == x) + 1
    super()[#pos]
    if institutions.keys().at(institutions.keys().position(y => y == x)) != authors.at(person).affiliations.last()  {super()[,]}
  }#if person != authors.keys().last(){","}]
]

v(0.25em)

// Print Corresponding Email
text(font: "Iosevka", size: 0.85em)[Corresponding Author: #corresponding_email]
v(1em)

// Print Institutions
for (i,inst) in institutions.keys().enumerate() {
  let pos = institutions.keys().position(y => y == inst) + 1
  text(0.8em, rgb("#222"))[#super()[#text(size: 1.2em, weight: "bold")[#pos ]]#institutions.at(inst)]
    if institutions.keys().last() != inst {text(size: 0.7em)[\ ]}
  }
}

#let showbibliography(file) = {
  pagebreak()
  columns(2)[
  #text(size: 0.8em)[
    #bibliography(file, title: "References", style: "springer-basic-author-date")
  ]]
}

#let note(author, body, color : orange) = {
  text(color)[(Note from #author: #body)]
}
