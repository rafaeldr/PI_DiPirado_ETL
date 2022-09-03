CREATE DEFINER=`root`@`localhost` PROCEDURE `VerificaInteracao`(IN nome1 VARCHAR(250), IN nome2 VARCHAR(250))
BEGIN
	
    DROP TEMPORARY TABLE IF EXISTS paTable1, paTable2;
  
	-- Init variables
	SET @id1 = 0;
	SET @id2 = 0;
	SET @origem1 = 0;
	SET @origem2 = 0;
	SET @count1 = 0;
	SET @count2 = 0;

    SELECT	id_principal, tipo_origem, count(id_principal) INTO @id1, @origem1, @count1
    FROM	bigtable_nomes 
    WHERE	nome LIKE CONCAT('%', nome1, '%')
    LIMIT 1;
	
    SELECT	id_principal, tipo_origem, count(id_principal) INTO @id2, @origem2, @count2
    FROM	bigtable_nomes 
    WHERE	nome LIKE CONCAT('%', nome2, '%')
    LIMIT 1;

    -- Recuperar IDs a partir da origem correta
	-- 1 KEGGDrug
	-- 2 anvisa_nome
	-- 3 anvisa_principioativo

    -- Origem 1: KEGGDrug
    IF @origem1 = 1 THEN
		CREATE TEMPORARY TABLE paTable1
			SELECT @id1 AS keggdrug_id, 1 as matchingValue;
	END IF;
    IF @origem2 = 1 THEN
		CREATE TEMPORARY TABLE paTable2
			SELECT @id2 AS keggdrug_id, 1 as matchingValue;
	END IF;    
    
    -- Origem 3: Anvisa_PrincipioAtivo
    IF @origem1 = 3 THEN
		CREATE TEMPORARY TABLE paTable1
			SELECT	keggdrug_id, matchingValue
			FROM	keggdrug_anvisa
			WHERE	id_pAtivo = @id1;
	END IF;
    IF @origem2 = 3 THEN
		CREATE TEMPORARY TABLE paTable2
			SELECT	keggdrug_id, matchingValue
			FROM	keggdrug_anvisa
			WHERE	id_pAtivo = @id2;
	END IF;
    
    -- Origem 2: Anvisa_Nome
    IF @origem1 = 2 THEN
		CREATE TEMPORARY TABLE paTable1 
			SELECT	k.keggdrug_id , k.matchingValue
			FROM	anvisa_nome_principioativo as prod
				INNER JOIN anvisa_principioativo as pa ON prod.idPrincipio = pa.id_pAtivo
				INNER JOIN keggdrug_anvisa as k ON k.id_pAtivo = pa.id_pAtivo
			WHERE	prod.idProduto = @id1;
	END IF;
    
	IF @origem2 = 2 THEN
		CREATE TEMPORARY TABLE paTable2
			SELECT	k.keggdrug_id, k.matchingValue
			FROM	anvisa_nome_principioativo as prod
				INNER JOIN anvisa_principioativo as pa ON prod.idPrincipio = pa.id_pAtivo
				INNER JOIN keggdrug_anvisa as k ON k.id_pAtivo = pa.id_pAtivo
			WHERE	prod.idProduto = @id2;
	END IF;

	SELECT	kn1.name, k.keggdrug_id1, kn2.name, k.keggdrug_id2, (p1.matchingValue + p2.matchingValue)/2 as matchingValue
	FROM	paTable1 as p1
		CROSS JOIN paTable2 as p2
		INNER JOIN keggdrug_interacao k 
			ON  k.keggdrug_id1 =
				CASE 
					WHEN p1.keggdrug_id < p2.keggdrug_id THEN p1.keggdrug_id
					WHEN p2.keggdrug_id < p1.keggdrug_id THEN p2.keggdrug_id
				END
			AND k.keggdrug_id2 = 
				CASE 
					WHEN p1.keggdrug_id > p2.keggdrug_id THEN p1.keggdrug_id
					WHEN p2.keggdrug_id > p1.keggdrug_id THEN p2.keggdrug_id
				END
		INNER JOIN keggdrug_nome kn1
			ON kn1.keggdrug_id = k.keggdrug_id1
		INNER JOIN keggdrug_nome kn2
			ON kn2.keggdrug_id = k.keggdrug_id2;


END