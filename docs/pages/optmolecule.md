---
title: Molecule
parent: Usage
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Molecule

Molecule block needs to be set up in every calculation. You need to specify the `"title" : "molecule"` line in the block to specify that you are defining the molecule. Below are all the availble options.

## Main Block Options

{: .note-title }
> path
>
> __Descriprion__: Path to the system you want to calculate in the .xyz format.<br>
> __Datatype__: string<br>
> __Default__: None

{: .note-title }
> basis
>
> __Descriprion__: No need to explain this. Just pick one of the avilable bases.<br>
> __Datatype__: string<br>
> __Default__: None

{: .note-title }
> charge
>
> __Descriprion__: Charge of the system.<br>
> __Datatype__: int<br>
> __Default__: 0

{: .note-title }
> multi
>
> __Descriprion__: Multiplicity of the system.<br>
> __Datatype__: int<br>
> __Default__: 1
