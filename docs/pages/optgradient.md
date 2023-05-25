---
title: Gradient
parent: Usage
layout: default
nav_order: 3
---
{% include mathjax.html %}

# Gradient

Specify `"title" : "engrad"` if you want to perform the gradient calculation. You need to specify the molecule block before this calculation. Below are the available options.

## Main Block Options

{: .note-title }
>  method
>
> __Description__: Specifies the method used as a json block. You can for example fill this with a Hartree-Fock block with some of these options.<br>
> __Datatype__: json<br>
> __Default__: None

{: .note-title }
> numerical
>
> __Description__: Specifies the usage of analytical or numerical gradient. If you specify analytical gradient for a method that is not compatible, an error will be printed.<br>
> __Datatype__: bool<br>
> __Default__: false

{: .note-title }
> increment
>
> __Description__: Only used if numerical gradient is enabled. Specifies the increment in coordinates in numerical derivatives.<br>
> __Datatype__: double<br>
> __Default__: 0.005

{: .note-title }
> nthread
>
> __Description__: Only used if numerical gradient is enabled. Specifies the number of cores the calculatio will use.<br>
> __Datatype__: int<br>
> __Default__: 1

{: .note-title }
> print
>
> __Description__: Some printing options.
>
> {: .note-title }
> > iter
> >
> > __Description__: Whether or not to print the iterations in numerical gradient calculation. Does nothing if analytical gradient is specified.<br>
> > __Datatype__: bool<br>
> > __Default__: true<br>
>
> {: .note-title }
> > method
> >
> > __Description__: Whether or not to print the used method results.<br>
> > __Datatype__: bool<br>
> > __Default__: true<br>
