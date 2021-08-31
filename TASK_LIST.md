# Sammanfattning av planen
- [x] Sätt upp en github,
- [x] Lägg till en README fil,
- [x] Sätt upp grenar eller "branches" på githuben,
- [x] Skapa en klassisk mappstruktur på Github,
- [x] Generera ett sfärisk nät via python mha gmsh,
- [x] Generera ett sfärisk nät med ett sfärisk hål via python mha gmsh,
- [x] Implementera så att meshen kan ses via paraview och python,
- [x] Fixa så att msh-näten som genererats av Gmsh kan konverteras till xdmf-nät som kan läsas in av Dolfin och Fenics vilket görs med hjälp av meshio,
- [ ] Skriv ett enkelt skript som läser in den konverterade xdmf-meshen via Fenics eller Dolfin,
- [ ] Sätt upp en variationsformulering,
- [ ] Implementera ett adaptivt tidssteg, 
- [ ] Reproducera fig. 4.3-4.6 i Chaplain (2001),
- [ ] Lägg till ett eller flera budscars. Beräkna perturberade egenvärden och undersök om de lämnar det instabila intervallet,
- [ ] Bekräfta de teoretiska resultaten genom att lösa systemet numeriskt,
- [ ] Inför en växande sfär med radie R(t),
- [ ] Undersök mönsterbildning på den växande sfären,
- [ ] Undersök hur närvaron av budscars påverkar möjligheten för mönsterbildning på den växande sfären,
- [x] Såhär kryssar man i en ruta. 



## Branches:
Här är en väldigt bra [länk](https://thenewstack.io/dont-mess-with-the-master-working-with-branches-in-git-and-github/) för hur man jobbar med grenar eller "branches" på github.

Jag har gjort en gren för varje deltagare i projektet:
1. Johannes's branch,
2. Philip's branch.
  
En viktig tumregel är att *allt arbete görs på de individuella grenarna* och sedan så slår man ihop arbetetet ("merge") vilket sammanfattas på main-branchen. Det finns två undantag till regeln om att inget arbete görs på main-grenen:

1. Filen README.md modifieras endast på main-grenen,
2. Filen TASK_LIST.md (alltså denna fil) modifieras likaså endast på main-grenen. 
