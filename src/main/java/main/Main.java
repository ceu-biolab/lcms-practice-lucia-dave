package main;

import lipid.LipidScoreUnit;
import org.drools.ruleunits.api.RuleUnitInstance;
import org.drools.ruleunits.api.RuleUnitProvider;

import java.util.List;

public class Main {

    public static void main(String[] args) {
        LipidScoreUnit lipidScoreUnit = new LipidScoreUnit();

        RuleUnitInstance<LipidScoreUnit> instance = RuleUnitProvider.get().createRuleUnitInstance(lipidScoreUnit);

        try {
            // TODO INTRODUCE THE CODE IF DESIRED TO INSERT FACTS AND TRIGGER RULES
            instance.fire();
            // TODO INTRODUCE THE QUERIES IF DESIRED


        } finally {
            instance.close();
        }
    }
}
