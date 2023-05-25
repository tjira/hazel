---
title: Molecular Dynamics
parent: Usage
layout: default
nav_order: 5
---
{% include mathjax.html %}

# Molecular Dynamics

Specify `"title" : "dynamics"` if you want to perform molecular dynamics. You need to specify the molecule block before this calculation. Below are the available options.

## Main Block Options

{: .note-title }
>  engrad
>
> __Description__: Specifies the gradient used as a json block. Pleas refer to the gradient options.<br>
> __Datatype__: json<br>
> __Default__: None

{: .note-title }
> timestep
>
> __Description__: Specifies a timestep for the dynamics.<br>
> __Datatype__: double<br>
> __Default__: None

{: .note-title }
> steps
>
> __Description__: Specifies number of steps.<br>
> __Datatype__: int<br>
> __Default__: None

{: .note-title }
> thermo
>
> __Description__: Enables a thermostat. If not specified, NVE dynamics is performed.
>
> {: .note-title }
> > title
> >
> > __Description__: Name of the thermostat. Currently avalable thermostat is berendsen.<br>
> > __Datatype__: string<br>
> > __Default__: None<br>
>
> {: .note-title }
> > temp
> >
> > __Description__: Temperature, that the thermostat will try to keep constant.<br>
> > __Datatype__: double<br>
> > __Default__: 298.15<br>
>
> {: .note-title }
> > q
> >
> > __Description__: Power of the temperature fraction. Only used in Berendsen thermost.<br>
> > __Datatype__: double<br>
> > __Default__: 0.5<br>
