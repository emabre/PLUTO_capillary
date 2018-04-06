# Avanzamento nel tempo delle equazioni in PLUTO

## Come è fatto l'avanzamento delle parti paraboliche
`parabolic_rhs.c / parabolicRHS()`: cacola il membro destro per algoritmi "split" (STS, RKC)
`parabolic_flux.c / parabolicFlux()` : calcola il "flusso" parabolico per algoritmo esplicito
Notare che sia `parabolicRHS()` sia  `parabolicFlux()` chiamano `ViscousFlux(..)`, `TC_Flux(..)`, `ResistiveFlux(..)` (ma ovviamente usano i risultati in modo diverso!)

## Come è fatto il ciclo principale di avanzamento delle equazioni
Dentro `main.c` si chiama `Integrate(..)`.
Integrate fa tutto il cuore dell'algoritmo:
A step alterni fa:
`AdvanceStep(..)`[*step non operator-split*] -> `split_source.c/SplitSource(..)`[*step operator-split*]
e
`split_source.c/SplitSource(..)`[*operator-split*] -> `AdvanceStep(..)`[*non operator-split*]

__Attenzione__:
(i) In caso non si stia usando l'agoritmo STS, nè RKC, nè cooling, `SplitSource(..)` non fa alcunchè!
(ii) Non confondere *operator-split* con *dimensional-splitting*: il primo si riferisce alla separazione degli operatori, il secondo (che è un caso particolare del primo) degli addendi dati dalle diverse dimensioni spaziali.

### Advance step
Essenzialemente, nel caso di metodi RKx,
`rk_step.c/AdvanceStep(..)`, quando avanza le equazioni usa `UpdateStage(..)` (contornato dalle opportune trasformazioni da variabili primitive a conservative e dalle pesature delle variabili secondo i coefficienti dati dal metodo di RK).
####`UpdateStage()`
Cosa fa:
+ Calcola $\vec J$ e $T$;
+ Chiama `States(..)`:
  - Calcola gli stati destro e sinistro alle interfacce, secondo un qualche algoritmo(ce ne sono vari disponibili))
+ Chiama *`Riemann(states,...)`* [puntatore ad un Riemann solver]:
  - Calcola il flusso iperbolico e lo salva in un campo della struttura `states` (penso `states.flux` e `states.press`)
+ **[se l'avanzamento dei termini parabolici è esplicito]** Chiama `ParabolicFlux(..)`
  - Chiama le tre funzioni per i flussi parabolici, che scrivono sui campi di `states` che sono i flussi parabolici di loro competenza:
    - `ViscousFlux(..)` -> `states.visc_flux`
    - `TC_Flux(..)` -> `states.tc_flux`
    - `ResistiveFlux(..)` -> `states.res_flux`
  - Aggiunge al flusso in `states` i tre flussi parabolici:
    - `states.flux` $+=$ `states.visc_flux`, `states.tc_flux`, `states.res_flux`
+ Chiama `RightHandSide(..)`:
  - costruisce il *right-hand-side* `state.rhs` per avanzare le equazioni, usando quantità che dipendono da varie cose, tra cui `states.flux` e `states.press`
+ Avanza nel tempo le variabili conservative facendo semplicemente:
  - `UU[..][..][..][..] += states.rhs[..][..]`

### SplitSource()
Considera radiative losses (qui li trascuro) e termini *parabolici*:
+ chiama `STS(..)` o `RKC(..)`
  - chiama `ParabolicRHS(..)` (la quale usa : `ViscousFlux(..)`, `TC_Flux(..)` e `ResistiveFlux(..)`)
  - Avanza alcune equazioni con Super-Time-Stepping a seconda di quali contributi parabolici sono definiti.
  `VISCOSITY` -> eq. momento, eq. energia
  `RESISTIVITY` ->  eq. induzione, eq. energia
  `THERMAL_CONDUCTION` -> eq. energia
