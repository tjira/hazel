---
title: Frequency Analysis
parent: Usage
layout: default
nav_order: 4
---
{% include mathjax.html %}

# Frequency Analysis

Specify `"title" : "freq"` if you want to perform the frequency analysis. You need to specify the molecule block before this calculation. Currently only numerical hessian is available. Below are the available options.

## Main Block Options

{: .note-title }
>  method
>
> __Description__: Specifies the method used as a json block. You can for example fill this with a Hartree-Fock block.<br>
> __Datatype__: json<br>
> __Default__: None

{: .note-title }
> increment
>
> __Description__: Specifies the increment in coordinates in numerical derivatives.<br>
> __Datatype__: double<br>
> __Default__: 0.005

{: .note-title }
> nthread
>
> __Description__: Specifies the number of cores the calculation will use.<br>
> __Datatype__: int<br>
> __Default__: 1

{: .note-title }
> scale
>
> __Description__: Scaling factor of the frequencies.<br>
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
> > __Description__: Whether or not to print the iterations in numerical hessian calculation.<br>
> > __Datatype__: bool<br>
> > __Default__: true<br>
>
> {: .note-title }
> > hessian
> >
> > __Description__: Whether or not to print the used method results.<br>
> > __Datatype__: bool<br>
> > __Default__: true<br>
