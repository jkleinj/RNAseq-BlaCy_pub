## functions to match gene names, IDs and indices

## function returning gene names of gene IDs
geneIDName = function(id, order.names = FALSE, join = FALSE) {
	name = as.data.frame(gene_names.unique$SYMBOL[match(id, gene_names.unique$ENSEMBL)]);
	## correct for NA names: replace clu_rnd.mat.rownames with ID number
	name[is.na(name)] = sapply(id[is.na(name)], function(x) { substr(x, 14, 18) });
	if (order.names == FALSE) {
		if (join == FALSE) {
			colnames(name) = c("gene_name");
			return(name);
		} else {
			id.name = as.data.frame(cbind(id, name));
			colnames(id.name) = c("gene_ID", "gene_name");
			return(id.name);
		}
	} else {
		if (join == FALSE) {
			colnames(name) = c("gene_name");
			return(name[order(name)]);
		} else {
			join = cbind(id, name);
			id.name = as.data.frame(join[order(join[, 2]), ]);
			colnames(id.name) = c("gene_ID", "gene_name");
			return(id.name);
		}
	}
}

## function returning gene IDs of gene names
geneNameID = function(name, order.IDs = FALSE, join = FALSE) {
	id = as.data.frame(gene_names.unique$ENSEMBL[match(name, gene_names.unique$SYMBOL)]);
	if (order.IDs == FALSE) {
		if (join == FALSE) {
			colnames(id) = c("gene_ID");
			return(id);
		} else {
			name.id = as.data.frame(cbind(name, id));
			colnames(name.id) = c("gene_name", "gene_ID");
			return(name.id);
		}
	} else {
		if (join == FALSE) {
			id.order = as.data.frame(id(order(id)));
			colnames(id.order) = c("gene_ID");
			return(id.order);
		} else {
			join = cbind(name, id);
			name.id = as.data.frame(join[order(join[, 2]), ]);
			colnames(name.id) = c("gene_name", "gene_ID");
			return(name.id);
		}
	}
}

## function returning indices of gene IDs
geneIDIdx = function(dat.id, query.id, order.idx = FALSE) {
    idx = which(dat.id %in% query.id);
    if (order.idx == FALSE) {
        return(idx);
    } else {
        return(idx[order(idx)]);
    }
}

## function returning indices of gene names 
geneNameIdx = function(dat.name, query.name, order.idx = FALSE) {
    idx = which(dat.name %in% query.name);
    if (order.idx == FALSE) {
        return(idx);
    } else {
        return(idx[order(idx)]);
    }
}

