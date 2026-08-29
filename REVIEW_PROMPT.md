# Prompt: przegląd prezentacji (wygląd + błędy)

Skopiuj poniższy tekst do Claude'a jako pojedynczą wiadomość.

---

Przejrzyj dokładnie moją prezentację Beamer i zaproponuj oraz wprowadź poprawki — zarówno wizualne, jak i merytoryczne.

## Kontekst projektu

- Repozytorium: `C:\Users\Szymon\Documents\github\MathBio Warsaw`
- Prezentacja: podkatalog `Latex/`. Skrypty pomocnicze: `Julia/`.
- Wykład: "How Does a Hydra Decide Where to Grow a Head?" — pattern formation, Turing instability, mechanochemistry. Slajdy są po angielsku.
- Silnik: **LuaLaTeX** (`fontspec`, font Arial). Buduj z katalogu `Latex/`: `latexmk -lualatex main.tex`.
- `Latex/main.tex` to punkt wejścia, ale **większość sekcji jest zakomentowana** — odkomentuj wszystkie `\input{sections/...}`, żeby zbudować pełną prezentację, i przywróć oryginalny stan `main.tex` na końcu.
- Kolejność sekcji: `A_hydra`, `B_reaction_diffusion`, `C_local_dynamics`, `D_turing_instability`, `GiererMeinhardt`, `E_simulations`, `F_mechanochemistry`. Razem ~159 slajdów, ~11 700 linii.
- `Latex/preamble.tex` — ciemny motyw (`HydraNavy` tło, akcenty `HydraGold`/`HydraCyan`/`HydraViolet`/`HydraMagenta`), własny `frametitle` z breadcrumbem.
- `Latex/macros.tex` — własne makra zdefiniowane przez `\NewDocumentCommand`, m.in. `\HydraSection`, `\HydraSubsection`, `\titlelayout`. Przeczytaj je zanim cokolwiek zmienisz.
- Figury w `Latex/figures/` (podfoldery: `hydra`, `diffusion`, `gierer_meinhardt`, `local_dynamics`, `turing`, `mechanochemistry`, `experiments`).
- To jest repozytorium git. Zrób commit stanu wyjściowego przed zmianami, a poprawki commituj partiami, żeby dało się je cofnąć pojedynczo.

## Krok 1 — obejrzyj prezentację

Zbuduj pełny PDF, a potem **obejrzyj go stronę po stronie** narzędziem Read (PDF czytaj partiami po ≤20 stron — przejrzyj cały dokument, nie próbkę). Równolegle czytaj źródła `.tex`, żeby wiedzieć, który kod odpowiada za który slajd.

Nie opieraj oceny wyglądu wyłącznie na kodzie — najpierw zobacz wyrenderowany slajd.

## Krok 2 — na co zwrócić uwagę

**Kompozycja i układ**
- treść wychodząca poza slajd, nachodzące na siebie elementy, ucięte wzory lub obrazki
- slajdy przepełnione (za dużo treści) i slajdy świecące pustkami
- niespójne marginesy, wyrównanie kolumn, wyśrodkowanie
- proporcje i rozdzielczość figur; czy podpisy są czytelne przy rozmiarze wyświetlania

**Typografia**
- W źródłach jest **~760 ręcznych `\fontsize{}{}`** i **~720 ręcznych `\vspace{}`** (często ujemnych). To główny dług projektu. Zidentyfikuj powtarzające się wartości i zaproponuj zestaw semantycznych makr (np. `\hydraCaption`, `\hydraFormula`, `\hydraLead`) w `macros.tex`, żeby rozmiary były spójne i zmieniane w jednym miejscu.
- za małe fonty (poniżej ~7 pt jest nieczytelne z sali), niespójne rozmiary dla tej samej roli tekstu
- rzeki, wdowy/sieroty, dziwne łamanie wierszy, nadużycie `\mbox`

**Kolor i kontrast**
- czy paleta jest stosowana konsekwentnie (czy dany kolor zawsze znaczy to samo)
- kontrast tekstu do tła `HydraNavy` — sprawdź zwłaszcza `HydraMuted` na małych rozmiarach
- czy kolory są rozróżnialne dla osób z daltonizmem; czy kolor nie jest jedynym nośnikiem informacji

**Spójność**
- format tytułów slajdów, styl breadcrumbów, jednolity układ slajdów sekcyjnych
- notacja matematyczna: te same symbole dla tych samych wielkości w całej prezentacji
- styl list, podpisów, wyróżnień, sposób podawania źródeł pod figurami

**Błędy**
- ostrzeżenia z `main.log`: Overfull/Underfull `\hbox`, undefined references, brakujące pliki, missing characters
- literówki i błędy językowe w angielskim tekście
- błędy w matematyce: indeksy, znaki, niespójne oznaczenia, złe odmiany operatorów
- błędy merytoryczne w treści biologicznej/matematycznej
- nieużywane pliki, martwy zakomentowany kod, zdublowane slajdy

## Krok 3 — najpierw raport, potem zmiany

**Zanim cokolwiek zmienisz**, wypisz mi **pełną, ponumerowaną listę wszystkich znalezionych problemów**. Dla każdego podaj:

1. numer strony PDF i plik z numerem linii (`Latex/sections/X.tex:123`)
2. na czym polega problem
3. konkretną proponowaną poprawkę
4. wagę: **krytyczny** (uszkodzony wygląd / błąd merytoryczny) / **istotny** (zauważalnie gorszy odbiór) / **kosmetyczny**

Pogrupuj listę tematycznie. Nie skracaj — chcę zobaczyć wszystko, także drobiazgi.

## Krok 4 — wprowadź poprawki

Po przedstawieniu listy wprowadź poprawki, w kolejności: krytyczne → istotne → kosmetyczne.

Zasady:
- **Nie zmieniaj treści merytorycznej** (twierdzeń, wzorów, wniosków naukowych) bez wyraźnego oznaczenia. Jeśli uważasz, że coś jest merytorycznie błędne — zgłoś to i zapytaj, zamiast poprawiać po cichu.
- Zachowaj istniejącą paletę i ogólny charakter motywu. Chodzi o dopracowanie, nie o przeprojektowanie od zera.
- Pracuj sekcja po sekcji, przebudowując PDF po każdej sekcji i sprawdzając wizualnie, że zmiany dały zamierzony efekt i niczego nie zepsuły.
- Poprawki systemowe (makra na rozmiary czcionek, ujednolicone odstępy) rób w `macros.tex`/`preamble.tex`, a nie przez powielanie łatek w każdej sekcji.
- Jeśli poprawka wymaga decyzji projektowej, na którą są dwie sensowne odpowiedzi — pokaż mi obie zamiast zgadywać.

Na końcu: przebuduj pełny PDF, potwierdź, że kompiluje się bez błędów, podsumuj co zostało zmienione i wypisz to, czego świadomie nie ruszyłeś (wraz z powodem).
