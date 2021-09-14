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
- [x] Skapa ett "launch"-skript där alla simuleringar körs ifrån,
- [x] Skapa ett "toolbox"-skript där alla FEM funktioner specificeras och som "launch"-skriptet baseras på,
- [x] Lägg till en funktion i toolboxen som beräknar kritiska värden på d och gamma (som kallas för d\_c och gamma\_c) samt steady state värdena i Schnackenberg modellen
- [x] Lägg till en funktion i toolboxen som läser in beräkningsnäten (dvs mesharna),
- [x] Lägg till en funktion i toolboxen som definierar alla funktionsrum i variationsformuleringen,
- [x] Lägg till en funktion i toolboxen som interpolerar ner steady staten på Hilbert rummet,
- [x] Lägg till en funktion som interpolerar ner en normalfördelad störning på Hilbert rummet,
- [x] Lägg till en funktion i toolboxen som sätter initialvillkoren baserat på en "dof-map" (dvs baserat på subdomänerna),
- [x] Fixa så att output filerna som sparas ner (dvs pvd och vtu) sparas ner i undermappar till mappen "Output" som är namngivna efter parametrarna i varje körning, 
- [x] Lägg till pvd och vtu filerna i .gitignore så att dessa inte pushas till github (dessa datafiler är enormt stora och sparas med fördel på Dropbox istället då github endast är till för att få tillgång till koden),
- [x] Spara ner initialvillkoren i en logiskt namngiven mapp i "Output" mappen och kontrollera dem i ParaView,
- [x] Skriv en README fil i Output mappen som beskriver hur filerna är namngivna,
- [x] Kontakta och svara en av grundarna till Gmsh vid namn Professor Christophe Geuzaine angående hur näten genereras,
- [x] Testa Christophes implementering,
- [x] Installera den nyaste versionen av Gmsh 4.8.4 via anaconda vilket innebär att skapa en ny conda omgivning vid namn "gmsh\_latest\_version" (detta innebär att vi kommer ha en miljö för FEniCS vid namn "fenicsproject" och en omgivning för Gmsh vid namn "gmsh\_latest\_version"),
- [x] Lägg till "Physical regions" i den nya meshen baserad på implementeringen som Christophe gjorde,
- [x] Lägg till färger i denna mesh och uppdatera vår figur av meshen,
- [x] Baserat på Christophes kod, skapa en mesh utan några hål,
- [x] Baserat på Christophes kod, skapa en mesh med ett hål,
- [x] Baserat på Christophes kod, skapa en mesh med två hål,
- [x] Baserat på Christophes kod, skapa en mesh med fem hål,
- [x] Uppdatera koverteringsskriptet så att näten med 0, 1, 2 och 5 hål respektive konverteras från msh format till xdmf format,
- [x] Importera xdmf näten med 0, 1, 2 och 5 hål i FEniCS och beräkna arean av alla ytor samt summan av ytorna. Kontrollera så att den totala arean blir ungefär 4*pi vilket är den analytiska ytarean av enhetssfären,
- [x] Kontakta och svara Anders Logg om FEniCS samt Carls kontaktuppgifter. Anders sa mycket riktigt att implementera PDEer på mångfalder med kurvatur kräver lite extra eftertanke och vidare sa han att den artikel som jag bifogade "[Automating the solution of PDEs on the sphere and other manifolds in FEniCS 1.2 ](https://gmd.copernicus.org/articles/6/2099/2013/)" av Marie Rognes från 2013 är perfekt att använda som källa. Det som är bra med denna artikel är att det finns massa exempelskript där vi antagligen kan låna kod för att ta hänsyn till den utåtriktade normalen på sfären, och Anders sa till och med att jag alltid kunde kontakta Marie om det behövdes (hon hade tydligen varit Anders gamla doktorand), 
- [ ] Skriv en funktion i toolboxen som sätter upp variationsformuleringen,
- [ ] Skriv en funktion i toolboxen som löser PDE problemet med FEM i rummet och FD i tiden,
- [ ] Implementera ett adaptivt tidssteg i den sistnämnda funktionen, 
- [ ] Reproducera fig. 4.3-4.6 i Chaplain (2001),
- [ ] Kör samma simuleringar som i Chaplain fast nätet som har ett hål,
- [ ] Kör samma simuleringar med parametrar som varierar spatiellt dvs med förstärkning av aktivering eller inaktivering i regionen intill hålet,
- [ ] Kör om Chaplains simuleringar fast med en växande sfär dvs där gamma varierar med tiden,
- [ ] Kör om simuleringarna med nätet med ett hål för en växande sfär,
- [ ] Kör om simuleringarna med nätet med ett hål för en växande sfär och med lokal förstärkning av parametrarna,
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
