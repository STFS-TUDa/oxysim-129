# Cloud Handler
↑ [OxySim-129](./../../README.md)

[CloudHandler.H](./CloudHandler.H)  
[CloudHandler.C](./CloudHandler.C)

## Overview 

The CloudHandler manages `carbonaceousCloud` instances at runtime. It handles cloud evolution, accumulates particle-to-gas exchange terms, and writes the clouds.

Multiple clouds of the same type can be combined. A dictionary with similar entries needs to be provided, where each key is the name of a cloud (free to choose, used e.g. in `deleteParticlesDict`, see below) and its value is the cloud type:
```
clouds
{ 
        biomassCloud    carbonaceousCloud;
        coalCloud       carbonaceousCloud;
}
```

Because clouds are only known by name and type at runtime (read from this dict, not hardcoded), the solver itself only ever calls a single accumulation function per exchange term — e.g. `cloudHandler.SrhoClouds()` in [`rhoEqn.H`](./../../applications/solvers/flameletCloudFoam/rhoEqn.H) or `cloudHandler.SUClouds(U, false)` in [`UEqn.H`](./../../applications/solvers/flameletCloudFoam/UEqn.H) — regardless of how many clouds are configured. Adding, removing, or renaming a cloud is a one-line change to the `clouds` dict; it never touches solver code or needs a recompile.

For an overview of the functions refer to [CloudHandler.H](./CloudHandler.H).

## Templating

Due to the structure of OpenFOAM with inheritance of templates, "specialized" functions are used to loop over the lists of clouds. These "specialized" functions are created with templates.

## Particle deletion

In large configurations, it may be useful to delete particles that have fully reacted or that stick to walls. This can reduce computation time significantly.
Inside the deleteParticlesDict, subdictionaries can be added for multiple deletion groups. Inside the deletion groups you can define further dictionaries to add or subtract particles (of a certain cloud) based on a criterion. The addition or subtraction are done in order.

Deletion occurs by default at every output time as defined by `writeInterval`. Set `deleteParticlesInterval <value>` to specify a custom deletion frequency.

In this example, particles are deleted every 0.001 seconds. First, all particles below the z coordinate of -0.599 are added to the deletion and then all particles whose temperature is higher than 700 K are removed from the deletion. This means particles below 700 K below -0.599 m are deleted.

```
deleteParticles true;

deleteParticlesInverval 0.001;

deleteParticlesDict
{
    deleteBottom
    {
    cloudName coalCloud;
    height
        {
            action  add;
            field   z;
            operator  less;
            value  -0.599;
        }
    temperature
        {
            action  subtract;
            field   T;
            operator  greater;
            value    700;
        }
    }
}
```

Currently available fields are:
- x, y, z
- rz (radial distance from z-axis)
- Uz (radial velocity around z-axis)
- Umag (velocity magnitude)
- T

New fields can be added in the code. Only particle variables that are accessible from the outside can be used.