//#import "authors.typ": *

#let preprint(
  metadata: "metadata.json",
  doc,
) = {
  set text(font: "Libertinus Serif", size: 12pt)
  show raw: set text(font: "Libertinus Sans")
  set par(spacing: 3em, justify: false)
  set page(paper: "us-letter", margin: 1.2in)

  show math.equation: set text(font: "Libertinus Math")


  set table.hline(stroke: .6pt)
  show figure.caption: it => {
    context {
      set align(left)
      set par(leading: 0.4em, hanging-indent: 0pt, justify: false)
      (
        text(10pt, it.supplement + " " + it.counter.display() + "\n", weight: "semibold", font: "Libertinus Sans")
          + text(9pt, it.body, font: "Libertinus Sans", luma(20%))
      )
    }
  }
  show figure.where(
    kind: table,
  ): set figure.caption(position: top)


  show heading.where(
    level: 1,
  ): it => block(width: 100%)[
    #v(1em)
    #set text(18pt, weight: 600, font: "Libertinus Sans", )
    // #smallcaps(it.body)
    #it.body
    #v(1em)
  ]

  show heading.where(
    level: 2,
  ): it => block(width: 100%)[
    #v(1em)
    #set text(14pt, weight: "semibold", fill: black.lighten(15%), font: "Libertinus Sans")
    #it.body
    #v(1.0em)
  ]

  show heading.where(
    level: 3,
  ): it => block(width: 100%)[
    #v(1em)
    #set text(fill: black.lighten(30%), weight: "semibold", font: "Libertinus Sans")

    #it.body
    #v(1.0em)
  ]

  let titlepage(data) = block[
    #v(2fr)
    #set par(spacing: 1em, leading: 0.3em, justify: false)
    #text(data.title, size: 30pt, font: "Libertinus Sans", weight: "regular")
    #v(1fr)

    
    #for (author) in data.authors {
      text(author.name, font: "Libertinus Sans")
  
      v(-0.4em)
      text(size: 0.75em, fill: black.lighten(40%), style: "italic", author.institution)
      if "email" in author {
        v(-0.6em)
        text(style: "italic", size: 0.7em, author.email)
      } 
      linebreak()
      linebreak()
    }
    #v(5fr)
  ]

  // Editing marks
  let add(body) = text(fill: rgb(0, 100, 0))[#body]
  let change(body) = text(fill: rgb(0, 100, 100))[#underline(body, stroke: rgb(0, 90, 90))]
  let cut(body) = text(fill: rgb(150, 150, 150))[#strike(body, stroke: rgb(100, 0, 0))]

  // let add(body) = text()[#body]
  // let change(body) = text()[#body]
  // let cut(body) = []

  show "TK": text(weight: "bold", font:"Libertinus Sans", fill: rgb("#e08619"))[TK]
  show "REF": text(weight: "bold", font:"Libertinus Sans", fill: rgb("#c6218c"))[REF]

  titlepage(json(metadata))
  pagebreak()

  
  set page(header: [
    #set text(font: "Libertinus Sans", size: 11pt, rgb("#333"))
    Catchen _et al._
    #h(1fr)
    _Projecting bumble bee pollination network disassembly_
  ])

  set page(footer: context [
    #set text(font: "Libertinus Sans", size: 10pt, rgb("#333"))
    Last update: #datetime.today().display()
    #h(1fr)
    Page
    #counter(page).display(
      "1 of 1",
      both: true,
    )
  ])

  [
    *Abstract*: Climate warming and land-use change are reshuffling the distribution of life on Earth. This change is altering the structure of species interaction networks, which ultimately enable the persistence of biodiversity and ecosystem services. Forecasting change in species interactions is a central challenge for biodiversity conservation, but there are numerous methodological challenges associated with spatiotemporally explicit mapping interactions because these interactions form networks that intrinsically vary in space and time. Here we project how interaction networks are rewired over time by integrating species distribution models with in-situ plant-pollinator interaction data to map expected change in a plant-pollinator network consisting of 13 bumble bee species and 157 plant species they forage from in the southern Rocky Mountains of Colorado, a region where strong elevation gradients could drive spatial mismatch under climate change. Models project the "vertical disassembly" of interaction networks, where elevational range shifts lead to increasingly large spatial mismatches under more extreme climate warming scenarios. Our models identify hotspots of change where up to 50% of the total number of interactions in the whole system are lost, often outpacing the arrival of new interactions. These results demonstrate the utility of species distribution projections in mapping the impact of global change on interaction networks.

    *Keywords*: pollination, interaction networks, biogeography, network rewiring, species distribution models, climate projections
  ]

  pagebreak()

  set par(leading: 18pt)
  set par.line(numbering: n => text(size: 10pt, font: "Libertinus Sans", luma(60%))[#n])
  set math.equation(numbering: "(1)")

  doc
}