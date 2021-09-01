# Sammanfattning av planen
- [x] Sätt upp en github,
- [x] Lägg till en README fil,
- [x] Sätt upp grenar eller "branches" på githuben,
- [x] Skapa en klassisk mappstruktur på Github,
- [x] Generera ett sfärisk nät via python mha gmsh,
- [x] Generera ett sfärisk nät med ett sfärisk hål via python mha gmsh,
- [x] Implementera så att meshen kan ses via paraview och python,
- [x] Fixa så att msh-näten som genererats av Gmsh kan konverteras till xdmf-nät som kan läsas in av Dolfin och Fenics vilket görs med hjälp av meshio,
- [x] Skriv ett enkelt skript som läser in den konverterade xdmf-meshen via Fenics eller Dolfin,
- [x] Kontrollera att rätt nät har lästs in genom att jämföra med enhetssfärens area och justera nogrannheten i nätet baserat på detta,
- [x] Skriv ett yml skript för att kunna köra projektet med alla viktiga paket i (MODIFERAD: Tog bort yml filen i slutändan ty det är lättare att modifera (mer specifikt lägga till paket) till den conda omgivning som fenicsprojektet har skapat),
- [x] Fixa så att installationen av alla paket (fenics, gmsh och paraview) fungerar enkelt med anaconda,
- [x] Skriv ner instruktioner i README filen om hur man installerar all paket med hjälp av anaconda,
- [x] Skapa en ny markdown fil, kalla den "VERSIONS\_OF\_PACKAGES.md", i main branchen som listar alla paket i omgivningen "fenicsproject",
- [x] Bestäm täthet i näten (resultat theta_step=pi/350) genom att jämföra den numeriska arean (A=12.5661) med den analytiska ytarean för enhetssfären som är A=4*pi=12.5663...
- [ ] I samma skript som i föregående punkt sätt domänspecifika initialvillkor som en perturbation kring steady state,
- [ ] Spara ner initialvillkoren i en logiskt namngiven mapp i "Output" mappen och studera dem i ParaView,
- [ ] Sätt upp en variationsformulering,
- [ ] Implementera ett adaptivt tidssteg, 
- [ ] Reproducera fig. 4.3-4.6 i Chaplain (2001),
- [ ] Kör simuleringar med nätet med ett hål,
- [ ] Kör simuleringar med parametrar som varierar spatiellt dvs med förstärkning av aktivering eller inaktivering i regionen intill hålet,
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
