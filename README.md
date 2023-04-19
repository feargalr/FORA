# FORA
Feargal's Over Representation Analysis

## Under Construction ##

## Installation
```{r example}
library(devtools)
install_github("feargalr/FORA")
library(FORA)
```


## Example Usage

```{r example}
gene_list_to_test <- c("ENSG00000157764", "ENSG00000137845", "ENSG00000141510", "ENSG00000175899", "ENSG00000012048")
vector_all_detected_genes <- c("ENSG00000157764", "ENSG00000137845", "ENSG00000141510", "ENSG00000175899", "ENSG00000012048",
              "ENSG00000168484", "ENSG00000116761", "ENSG00000082898", "ENSG00000170781", "ENSG00000124145",
              "ENSG00000186106", "ENSG00000103169", "ENSG00000120688", "ENSG00000198947", "ENSG00000111422",
              "ENSG00000101266", "ENSG00000080824", "ENSG00000134250", "ENSG00000134222", "ENSG00000105974",
              "ENSG00000105819", "ENSG00000164118", "ENSG00000169189", "ENSG00000126945", "ENSG00000123609",
              "ENSG00000099991", "ENSG00000198931", "ENSG00000137942", "ENSG00000143602", "ENSG00000171791",
              "ENSG00000119918", "ENSG00000131469", "ENSG00000135046", "ENSG00000198804", "ENSG00000113810",
              "ENSG00000125488", "ENSG00000111679", "ENSG00000038295", "ENSG00000142156", "ENSG00000106495",
              "ENSG00000137869", "ENSG00000166503", "ENSG00000118699", "ENSG00000116001", "ENSG00000140538",
              "ENSG00000174574", "ENSG00000167133", "ENSG00000134108", "ENSG00000141232", "ENSG00000148175",
              "ENSG00000144455", "ENSG00000157823", "ENSG00000152683", "ENSG00000135297", "ENSG00000118685",
              "ENSG00000176495", "ENSG00000118729", "ENSG00000185315", "ENSG00000126003", "ENSG00000162552",
              "ENSG00000143793", "ENSG00000135446", "ENSG00000135756", "ENSG00000164005", "ENSG00000105198",
              "ENSG00000167478", "ENSG00000172316", "ENSG00000111201", "ENSG00000140442", "ENSG00000151010",
              "ENSG00000109321", "ENSG00000125968", "ENSG00000143549", "ENSG00000169814", "ENSG00000100234",
              "ENSG00000160237", "ENSG00000198286", "ENSG00000120509", "ENSG00000174047", "ENSG00000100678",
              "ENSG00000172378", "ENSG00000112851", "ENSG00000130600", "ENSG00000123668", "ENSG00000159266",
              "ENSG00000143977", "ENSG00000108468", "ENSG00000154122", "ENSG00000105329")

FORA(gene_list_to_test,vector_all_detected_genes)
```
