
## bp_foo  = with prifix and --fast
## bp_foo2 = with prefix but no --fast

SHOW TABLES LIKE "%attributelist%";

DROP TABLE bp_foo_attribute      ;
DROP TABLE bp_foo_attributelist  ;
DROP TABLE bp_foo_feature        ;
DROP TABLE bp_foo_interval_stats ;
DROP TABLE bp_foo_locationlist   ;
DROP TABLE bp_foo_meta           ;
DROP TABLE bp_foo_name           ;
DROP TABLE bp_foo_parent2child   ;
DROP TABLE bp_foo_sequence       ;
DROP TABLE bp_foo_typelist       ;




SELECT COUNT(*) AS bp_foo_attribute      FROM bp_foo_attribute      ;
SELECT COUNT(*) AS bp_foo_attributelist  FROM bp_foo_attributelist  ;
SELECT COUNT(*) AS bp_foo_feature        FROM bp_foo_feature        ;
SELECT COUNT(*) AS bp_foo_interval_stats FROM bp_foo_interval_stats ;
SELECT COUNT(*) AS bp_foo_locationlist   FROM bp_foo_locationlist   ;
SELECT COUNT(*) AS bp_foo_meta           FROM bp_foo_meta           ;
SELECT COUNT(*) AS bp_foo_name           FROM bp_foo_name           ;
SELECT COUNT(*) AS bp_foo_parent2child   FROM bp_foo_parent2child   ;
SELECT COUNT(*) AS bp_foo_sequence       FROM bp_foo_sequence       ;
SELECT COUNT(*) AS bp_foo_typelist       FROM bp_foo_typelist       ;




SHOW TABLES LIKE "bp_foo%";

SELECT COUNT(*) AS bp_foo_attribute      FROM bp_foo_attribute      ;
SELECT COUNT(*) AS bp_foo_attributelist  FROM bp_foo_attributelist  ;
SELECT COUNT(*) AS bp_foo_feature        FROM bp_foo_feature        ;
SELECT COUNT(*) AS bp_foo_interval_stats FROM bp_foo_interval_stats ;
SELECT COUNT(*) AS bp_foo_locationlist   FROM bp_foo_locationlist   ;
SELECT COUNT(*) AS bp_foo_meta           FROM bp_foo_meta           ;
SELECT COUNT(*) AS bp_foo_name           FROM bp_foo_name           ;
SELECT COUNT(*) AS bp_foo_parent2child   FROM bp_foo_parent2child   ;
SELECT COUNT(*) AS bp_foo_sequence       FROM bp_foo_sequence       ;
SELECT COUNT(*) AS bp_foo_typelist       FROM bp_foo_typelist       ;

SELECT COUNT(*) AS bp_foo2_attribute      FROM bp_foo2_attribute      ;
SELECT COUNT(*) AS bp_foo2_attributelist  FROM bp_foo2_attributelist  ;
SELECT COUNT(*) AS bp_foo2_feature        FROM bp_foo2_feature        ;
SELECT COUNT(*) AS bp_foo2_interval_stats FROM bp_foo2_interval_stats ;
SELECT COUNT(*) AS bp_foo2_locationlist   FROM bp_foo2_locationlist   ;
SELECT COUNT(*) AS bp_foo2_meta           FROM bp_foo2_meta           ;
SELECT COUNT(*) AS bp_foo2_name           FROM bp_foo2_name           ;
SELECT COUNT(*) AS bp_foo2_parent2child   FROM bp_foo2_parent2child   ;
SELECT COUNT(*) AS bp_foo2_sequence       FROM bp_foo2_sequence       ;
SELECT COUNT(*) AS bp_foo2_typelist       FROM bp_foo2_typelist       ;


## Problem tables...
# bp_foo2_attribute
# bp_foo2_feature
# bp_foo2_interval_stats


