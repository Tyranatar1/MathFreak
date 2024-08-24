// Normal Distribution Solver
function solveNormal() {
    const mu = parseFloat(document.getElementById('normal-mu').value);
    const sigma = parseFloat(document.getElementById('normal-sigma').value);
    const x = parseFloat(document.getElementById('normal-x').value);

    const pdf = (x, mu, sigma) => {
        return (1 / (sigma * Math.sqrt(2 * Math.PI))) * Math.exp(-Math.pow(x - mu, 2) / (2 * Math.pow(sigma, 2)));
    };

    const result = pdf(x, mu, sigma);
    document.getElementById('normal-result').innerText = `PDF Value: ${result}`;
}

// Binomial Distribution Solver
function solveBinomial() {
    const n = parseInt(document.getElementById('binomial-n').value);
    const p = parseFloat(document.getElementById('binomial-p').value);
    const k = parseInt(document.getElementById('binomial-k').value);

    const binomialCoeff = (n, k) => {
        let res = 1;
        for (let i = 0; i < k; i++) {
            res *= (n - i);
            res /= (i + 1);
        }
        return res;
    };

    const pmf = (n, k, p) => {
        return binomialCoeff(n, k) * Math.pow(p, k) * Math.pow(1 - p, n - k);
    };

    const result = pmf(n, k, p);
    document.getElementById('binomial-result').innerText = `PMF Value: ${result}`;
}

// Poisson Distribution Solver
function solvePoisson() {
    const lambda = parseFloat(document.getElementById('poisson-lambda').value);
    const k = parseInt(document.getElementById('poisson-k').value);

    const pmf = (lambda, k) => {
        return (Math.pow(lambda, k) * Math.exp(-lambda)) / factorial(k);
    };

    const factorial = (n) => {
        if (n === 0) return 1;
        let result = 1;
        for (let i = 1; i <= n; i++) {
            result *= i;
        }
        return result;
    };

    const result = pmf(lambda, k);
    document.getElementById('poisson-result').innerText = `PMF Value: ${result}`;
}

// Exponential Distribution Solver
function solveExponential() {
    const lambda = parseFloat(document.getElementById('exponential-lambda').value);
    const x = parseFloat(document.getElementById('exponential-x').value);

    const pdf = (x, lambda) => {
        return lambda * Math.exp(-lambda * x);
    };

    const result = pdf(x, lambda);
    document.getElementById('exponential-result').innerText = `PDF Value: ${result}`;
}

// Uniform Distribution Solver
function solveUniform() {
    const a = parseFloat(document.getElementById('uniform-a').value);
    const b = parseFloat(document.getElementById('uniform-b').value);
    const x = parseFloat(document.getElementById('uniform-x').value);

    const pdf = (x, a, b) => {
        return (x >= a && x <= b) ? 1 / (b - a) : 0;
    };

    const cdf = (x, a, b) => {
        if (x < a) return 0;
        if (x > b) return 1;
        return (x - a) / (b - a);
    };

    const resultPdf = pdf(x, a, b);
    const resultCdf = cdf(x, a, b);
    document.getElementById('uniform-result').innerText = `PDF Value: ${resultPdf}, CDF Value: ${resultCdf}`;
}

// Bernoulli Distribution Solver
function solveBernoulli() {
    const p = parseFloat(document.getElementById('bernoulli-p').value);
    const k = parseInt(document.getElementById('bernoulli-k').value);

    const pmf = (k, p) => {
        return k === 1 ? p : 1 - p;
    };

    const result = pmf(k, p);
    document.getElementById('bernoulli-result').innerText = `PMF Value: ${result}`;
}

// Gamma Distribution Solver
function solveGamma() {
    const k = parseFloat(document.getElementById('gamma-k').value);
    const theta = parseFloat(document.getElementById('gamma-theta').value);
    const x = parseFloat(document.getElementById('gamma-x').value);

    const gammaFunc = (z) => {
        // Using Gamma function approximation
        return Math.exp(z * Math.log(z) - z - Math.log(Math.sqrt(2 * Math.PI * z)));
    };

    const pdf = (x, k, theta) => {
        return (Math.pow(x, k - 1) * Math.exp(-x / theta)) / (Math.pow(theta, k) * gammaFunc(k));
    };

    const result = pdf(x, k, theta);
    document.getElementById('gamma-result').innerText = `PDF Value: ${result}`;
}

// Beta Distribution Solver
function solveBeta() {
    const alpha = parseFloat(document.getElementById('beta-alpha').value);
    const beta = parseFloat(document.getElementById('beta-beta').value);
    const x = parseFloat(document.getElementById('beta-x').value);

    const betaFunc = (a, b) => {
        return (gammaFunc(a) * gammaFunc(b)) / gammaFunc(a + b);
    };

    const gammaFunc = (z) => {
        // Using Gamma function approximation
        return Math.exp(z * Math.log(z) - z - Math.log(Math.sqrt(2 * Math.PI * z)));
    };

    const pdf = (x, alpha, beta) => {
        return (Math.pow(x, alpha - 1) * Math.pow(1 - x, beta - 1)) / betaFunc(alpha, beta);
    };

    const result = pdf(x, alpha, beta);
    document.getElementById('beta-result').innerText = `PDF Value: ${result}`;
}

// Negative Binomial Distribution Solver
function solveNegativeBinomial() {
    const r = parseInt(document.getElementById('negative-binomial-r').value);
    const p = parseFloat(document.getElementById('negative-binomial-p').value);
    const k = parseInt(document.getElementById('negative-binomial-k').value);

    const binomialCoeff = (n, k) => {
        let res = 1;
        for (let i = 0; i < k; i++) {
            res *= (n - i);
            res /= (i + 1);
        }
        return res;
    };

    const pmf = (r, p, k) => {
        return binomialCoeff(k + r - 1, k) * Math.pow(1 - p, k) * Math.pow(p, r);
    };

    const result = pmf(r, p, k);
    document.getElementById('negative-binomial-result').innerText = `PMF Value: ${result}`;
}
