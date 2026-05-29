# NEXT_STEPS.md — Co musisz zrobić ręcznie

Poniżej lista kroków, które wymagają Twojej ręcznej interwencji po przejęciu repozytorium.

---

## 1. Przejęcie repozytorium na konto dawidx1233

Repozytorium zostało utworzone na koncie `lenax04` (konto powiązane z tokenem GitHub w sandboxie). Aby przenieść je na konto `dawidx1233`:

1. Zaloguj się na GitHub jako `lenax04`.
2. Przejdź do: `https://github.com/lenax04/microbiome-pipeline/settings`
3. Przewiń do sekcji **Danger Zone** → **Transfer ownership**.
4. Wpisz `dawidx1233` jako nowego właściciela i potwierdź.
5. Zaktualizuj linki w `README.md`, `CITATION.cff` i `.github/workflows/ci.yml` z `lenax04` na `dawidx1233`.

---

## 2. Uzupełnienie danych osobowych

W następujących plikach wstaw swoje imię, afiliację i ORCID:

**`paper/paper.md`** (linia 3-4):
```
**Author:** Dawid [TWOJE NAZWISKO]
**Affiliation:** [1] [TWOJA INSTYTUCJA], [MIASTO], [KRAJ]
```

**`CITATION.cff`**:
```yaml
authors:
  - family-names: "[TWOJE NAZWISKO]"
    given-names: "Dawid"
    orcid: "https://orcid.org/[TWÓJ ORCID]"
```

Jeśli nie masz ORCID, zarejestruj się bezpłatnie na: https://orcid.org/register

---

## 3. Dodanie CI workflow przez GitHub UI

Token GitHub w sandboxie nie ma uprawnień do tworzenia plików workflow. Aby dodać CI:

1. Przejdź do: `https://github.com/dawidx1233/microbiome-pipeline/new/main`
2. Utwórz plik: `.github/workflows/ci.yml`
3. Wklej zawartość z pliku lokalnego: `/home/ubuntu/microbiome-pipeline/.github/workflows/ci.yml`

---

## 4. Rejestracja DOI przez Zenodo

Aby uzyskać trwały identyfikator DOI dla projektu:

1. Zaloguj się na https://zenodo.org (możesz użyć konta GitHub).
2. Przejdź do: **GitHub** → **Settings** → **Webhooks** → włącz integrację z Zenodo.
3. Lub ręcznie: **New Upload** → prześlij archiwum ZIP repozytorium.
4. Po uzyskaniu DOI zaktualizuj:
   - `CITATION.cff`: pole `doi`
   - `README.md`: badge DOI
   - `paper/paper.md`: sekcja References [1]

---

## 5. Zgłoszenie do czasopisma

Artykuł (`paper/paper.md` i `paper/paper.pdf`) jest sformatowany zgodnie z wymaganiami **GigaScience** (Oxford Academic).

Kroki do zgłoszenia:
1. Przejdź do: https://academic.oup.com/gigascience/pages/general_instructions
2. Sprawdź aktualne wymagania formatowania (Author Guidelines).
3. Uzupełnij brakujące dane: imię/afiliacja, ORCID, podziękowania (Acknowledgements), deklaracja konfliktu interesów.
4. Prześlij `paper/paper.pdf` przez system Editorial Manager.

Alternatywnie możesz zgłosić do **BMC Bioinformatics**: https://bmcbioinformatics.biomedcentral.com/submission-guidelines

---

## 6. Deponowanie danych (GigaScience / GigaDB)

GigaScience wymaga deponowania danych w **GigaDB** (http://gigadb.org):

1. Zarejestruj konto na GigaDB.
2. Prześlij pliki wynikowe z katalogu `results/` oraz dane testowe z `tests/data/`.
3. Uzyskaj DOI dla datasetu z GigaDB.
4. Dodaj sekcję **Data Availability** w `paper/paper.md` z linkiem do GigaDB.

---

## 7. Uruchomienie na serwerze zdalnym (partnerbdo.pl)

Klucz SSH dostarczony w zadaniu był uszkodzony (pola długości w strukturze binarnej były błędne). Aby uruchomić pipeline na serwerze:

1. Wygeneruj nowy klucz SSH: `ssh-keygen -t ed25519 -f ~/.ssh/microbiome_key`
2. Dodaj klucz publiczny do serwera: `ssh-copy-id -i ~/.ssh/microbiome_key.pub root@46.101.182.180`
3. Sklonuj repozytorium na serwerze: `git clone https://github.com/dawidx1233/microbiome-pipeline.git`
4. Zainstaluj Miniconda i Snakemake.
5. Uruchom pipeline: `snakemake --use-conda --cores 8`

---

## 8. Opcjonalne ulepszenia (PhD/stażowe portfolio)

Poniższe ulepszenia znacznie zwiększą atrakcyjność projektu:

- Dodanie wsparcia dla **long-read sequencing** (Oxford Nanopore / PacBio).
- Integracja z **QIIME2** jako alternatywny moduł amplikon.
- Moduł **machine learning** do biomarker discovery (Random Forest, XGBoost).
- Benchmark na publicznych danych HMP (Human Microbiome Project): https://hmpdacc.org/
- Benchmark na danych ENA/SRA: np. projekt PRJEB11419 (IBD microbiome).
- Dodanie **Nextflow** jako alternatywnego workflow managera.
